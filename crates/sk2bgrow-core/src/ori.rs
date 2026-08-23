//! Replication-origin handling: annotation lookup and a coarse grid search.
//!
//! The report's §6.1 argument is that TGT anchors carry real coordinates, so the
//! ori can be a *parameter* rather than something the sorted-regression trick
//! throws away. Two paths follow from that:
//!
//! * **ori known** — a DoriC / Ori-Finder annotation, or a `dnaA` position, is
//!   read from a table and used directly.
//! * **ori unknown** — the ori coordinate is searched jointly with the slope.
//!
//! This module owns the *coarse* search: a fast circular grid scan that returns
//! a starting point and a confidence. The refined piecewise V-shape MLE (with
//! the multi-fork plateau term) lives in `python/sk2bgrow/fit.py`, which is where
//! the report puts the statistics layer.

use std::collections::HashMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::anchor_db::AnchorDb;
use crate::error::{Result, Sk2bError};

/// Circular distance from `x` to `ori` on a genome of length `g`, in bp.
/// Ranges over `0 ..= g/2`; the maximum is reached at the terminus.
#[inline]
pub fn circular_distance(x: u64, ori: u64, g: u64) -> f64 {
    if g == 0 {
        return 0.0;
    }
    let d = (x as i64 - ori as i64).rem_euclid(g as i64) as u64;
    d.min(g - d) as f64
}

/// An optional prior on the ori position, e.g. from a `dnaA` annotation.
/// GRiD and SMEG both use `dnaA` as a QC anchor (report §6.1).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct OriPrior {
    pub position: u64,
    /// Standard deviation in bp. A curated oriC gets a tight prior; a `dnaA`
    /// position gets a loose one, since `dnaA` is near but not at the origin.
    pub sd: f64,
}

/// Result of a grid search.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OriFit {
    pub ori: u64,
    /// Terminus, i.e. the antipode of `ori`.
    pub ter: u64,
    /// Slope of log2(count) against circular distance; negative for a real
    /// replication gradient.
    pub slope: f64,
    pub intercept: f64,
    /// `-slope * genome_len / 2`: the log2 peak-to-trough ratio.
    pub log2_ptr: f64,
    pub r2: f64,
    pub sse: f64,
    /// Circular mean resultant length of the posterior over ori, in `[0, 1]`.
    /// 1.0 means the profile picks one position sharply; near 0 means the data
    /// are consistent with the ori being almost anywhere — the identifiability
    /// failure the report flags for slow growers (risk R4).
    pub confidence: f64,
    pub n_windows: usize,
}

impl OriFit {
    /// `PTR` on the linear scale.
    pub fn ptr(&self) -> f64 {
        self.log2_ptr.exp2()
    }
    /// A fit is only interpretable when the gradient runs downhill from the ori.
    /// An uphill slope means the search latched onto the terminus, or there is
    /// no replication signal at all.
    pub fn is_oriented_correctly(&self) -> bool {
        self.slope < 0.0
    }
}

/// Grid-search the ori on a circular genome.
///
/// `positions` are global coordinates of window midpoints and `values` the
/// matching log2 mean counts. `n_grid` candidate origins are tested evenly
/// around the circle; at each, log2(count) is regressed on circular distance and
/// the residual sum of squares recorded.
///
/// Returns `None` when there are fewer than three usable windows — two points
/// fit any line exactly and would report a meaningless R².
pub fn grid_search(
    positions: &[u64],
    values: &[f64],
    genome_len: u64,
    n_grid: usize,
    prior: Option<OriPrior>,
) -> Option<OriFit> {
    assert_eq!(
        positions.len(),
        values.len(),
        "positions and values must be parallel"
    );
    let usable: Vec<(u64, f64)> = positions
        .iter()
        .zip(values)
        .filter(|(_, v)| v.is_finite())
        .map(|(&p, &v)| (p, v))
        .collect();
    if usable.len() < 3 || genome_len == 0 {
        return None;
    }
    let n_grid = n_grid.max(4);
    let step = (genome_len as f64 / n_grid as f64).max(1.0);

    let mut best: Option<(f64, OriFit)> = None;
    let mut objective: Vec<(u64, f64)> = Vec::with_capacity(n_grid);

    for k in 0..n_grid {
        let ori = ((k as f64 * step) as u64) % genome_len;
        let (slope, intercept, sse, r2) = fit_line(&usable, ori, genome_len);
        // Only a downhill gradient is a replication profile; an uphill fit is
        // the same line read from the terminus and must not win.
        if !slope.is_finite() || slope >= 0.0 {
            objective.push((ori, f64::INFINITY));
            continue;
        }
        let penalty = match prior {
            Some(p) if p.sd > 0.0 => {
                let d = circular_distance(ori, p.position, genome_len) / p.sd;
                // Scaled to the residual variance so the prior competes on the
                // same footing as the data rather than swamping it.
                0.5 * d * d * (sse / usable.len() as f64)
            }
            _ => 0.0,
        };
        let obj = sse + penalty;
        objective.push((ori, obj));
        let fit = OriFit {
            ori,
            ter: (ori + genome_len / 2) % genome_len,
            slope,
            intercept,
            log2_ptr: -slope * (genome_len as f64 / 2.0),
            r2,
            sse,
            confidence: 0.0,
            n_windows: usable.len(),
        };
        // `Option::is_none_or` is stable only since 1.82; the workspace
        // MSRV is 1.74.
        #[allow(clippy::unnecessary_map_or)]
        if best.as_ref().map_or(true, |(b, _)| obj < *b) {
            best = Some((obj, fit));
        }
    }

    let (_, mut fit) = best?;
    fit.confidence = posterior_concentration(&objective, genome_len, usable.len());
    Some(fit)
}

/// Weighted least squares of `value ~ intercept + slope * circular_distance`.
fn fit_line(points: &[(u64, f64)], ori: u64, g: u64) -> (f64, f64, f64, f64) {
    let n = points.len() as f64;
    let xs: Vec<f64> = points
        .iter()
        .map(|&(p, _)| circular_distance(p, ori, g))
        .collect();
    let mx = xs.iter().sum::<f64>() / n;
    let my = points.iter().map(|&(_, v)| v).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (x, &(_, y)) in xs.iter().zip(points) {
        sxx += (x - mx) * (x - mx);
        sxy += (x - mx) * (y - my);
    }
    if sxx <= 0.0 {
        return (f64::NAN, my, f64::INFINITY, 0.0);
    }
    let slope = sxy / sxx;
    let intercept = my - slope * mx;
    let mut sse = 0.0;
    let mut sst = 0.0;
    for (x, &(_, y)) in xs.iter().zip(points) {
        let r = y - (intercept + slope * x);
        sse += r * r;
        sst += (y - my) * (y - my);
    }
    let r2 = if sst > 0.0 { 1.0 - sse / sst } else { 0.0 };
    (slope, intercept, sse, r2)
}

/// Turn the SSE profile into a posterior over ori and report its circular mean
/// resultant length.
///
/// Under a Gaussian error model the relative likelihood of two origins is
/// `exp(-n/2 · Δ log SSE)`, so the profile is converted with that transform
/// rather than an arbitrary temperature.
fn posterior_concentration(objective: &[(u64, f64)], g: u64, n_obs: usize) -> f64 {
    let min = objective
        .iter()
        .map(|&(_, o)| o)
        .fold(f64::INFINITY, f64::min);
    if !min.is_finite() || min <= 0.0 {
        return 0.0;
    }
    let mut sum_w = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for &(ori, obj) in objective {
        if !obj.is_finite() {
            continue;
        }
        let w = (-(n_obs as f64) / 2.0 * (obj / min).ln()).exp();
        if !w.is_finite() {
            continue;
        }
        let theta = 2.0 * std::f64::consts::PI * (ori as f64 / g as f64);
        sx += w * theta.cos();
        sy += w * theta.sin();
        sum_w += w;
    }
    if sum_w <= 0.0 {
        return 0.0;
    }
    ((sx / sum_w).powi(2) + (sy / sum_w).powi(2))
        .sqrt()
        .clamp(0.0, 1.0)
}

/// One row of an ori annotation table.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OriAnnotation {
    pub genome: String,
    pub position: u64,
    pub confidence: f64,
    pub source: String,
}

/// Read an ori annotation TSV.
///
/// Columns: `genome<TAB>position[<TAB>confidence][<TAB>source]`. Lines starting
/// with `#` are comments. This is the shape a DoriC or Ori-Finder export is
/// reduced to; keeping it a plain table means a curated origin, a `dnaA`
/// coordinate and a previous grid-search result are all loadable the same way.
pub fn read_annotations(path: &Path) -> Result<Vec<OriAnnotation>> {
    let text = std::fs::read_to_string(path).map_err(|e| Sk2bError::Io {
        path: path.into(),
        source: e,
    })?;
    let mut out = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 {
            return Err(Sk2bError::Format {
                path: path.to_path_buf(),
                msg: format!("line {}: expected at least genome<TAB>position", lineno + 1),
            });
        }
        let position: u64 = f[1].trim().parse().map_err(|_| Sk2bError::Format {
            path: path.to_path_buf(),
            msg: format!("line {}: '{}' is not a coordinate", lineno + 1, f[1]),
        })?;
        out.push(OriAnnotation {
            genome: f[0].trim().to_string(),
            position,
            confidence: f.get(2).and_then(|s| s.trim().parse().ok()).unwrap_or(1.0),
            source: f
                .get(3)
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| "annotation".into()),
        });
    }
    Ok(out)
}

/// Attach annotations to a database by genome name. Returns how many genomes
/// were matched; unmatched genomes keep `ori = None` and fall back to the grid
/// search.
pub fn attach(db: &mut AnchorDb, annotations: &[OriAnnotation]) -> usize {
    let map: HashMap<&str, &OriAnnotation> =
        annotations.iter().map(|a| (a.genome.as_str(), a)).collect();
    let mut hit = 0;
    for g in db.genomes.iter_mut() {
        if let Some(a) = map.get(g.name.as_str()) {
            if a.position >= g.genome_len {
                continue; // out of range for this assembly; ignore rather than wrap
            }
            g.ori = Some(a.position);
            g.ori_confidence = a.confidence;
            hit += 1;
        }
    }
    hit
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Windows on a circular genome with a clean V profile around `ori`.
    fn v_profile(g: u64, ori: u64, log2_ptr: f64, n: usize) -> (Vec<u64>, Vec<f64>) {
        let slope = -log2_ptr / (g as f64 / 2.0);
        let positions: Vec<u64> = (0..n).map(|i| (i as u64 * g) / n as u64).collect();
        let values = positions
            .iter()
            .map(|&p| 5.0 + slope * circular_distance(p, ori, g))
            .collect();
        (positions, values)
    }

    #[test]
    fn circular_distance_wraps() {
        let g = 1_000u64;
        assert_eq!(circular_distance(10, 990, g), 20.0);
        assert_eq!(circular_distance(990, 10, g), 20.0);
        assert_eq!(circular_distance(500, 0, g), 500.0, "terminus is g/2 away");
        assert_eq!(circular_distance(0, 0, g), 0.0);
    }

    #[test]
    fn grid_search_recovers_a_planted_ori() {
        let g = 4_600_000u64;
        let ori = 3_923_883u64; // E. coli K-12 oriC, as used in the report's simulation
        let (pos, val) = v_profile(g, ori, 1.0, 200);
        let fit = grid_search(&pos, &val, g, 256, None).unwrap();
        let err = circular_distance(fit.ori, ori, g);
        assert!(err < 2.0 * (g as f64 / 256.0), "ori off by {err} bp");
        assert!(
            (fit.log2_ptr - 1.0).abs() < 0.05,
            "log2 PTR was {}",
            fit.log2_ptr
        );
        assert!(fit.r2 > 0.99);
        assert!(fit.is_oriented_correctly());
        assert!(
            fit.confidence > 0.9,
            "a noiseless profile should be sharp, got {}",
            fit.confidence
        );
    }

    #[test]
    fn ptr_is_the_exponentiated_slope() {
        let g = 1_000_000u64;
        let (pos, val) = v_profile(g, 0, 2.0, 100);
        let fit = grid_search(&pos, &val, g, 128, None).unwrap();
        assert!((fit.ptr() - 4.0).abs() < 0.1, "PTR was {}", fit.ptr());
        assert_eq!(fit.ter, g / 2);
    }

    #[test]
    fn a_flat_profile_reports_low_confidence() {
        let g = 1_000_000u64;
        let pos: Vec<u64> = (0..100).map(|i| i * g / 100).collect();
        // Deterministic pseudo-noise around a flat mean: no gradient to find.
        let val: Vec<f64> = (0..100)
            .map(|i| 5.0 + ((i * 7919 % 13) as f64 - 6.0) * 0.01)
            .collect();
        let fit = grid_search(&pos, &val, g, 128, None).unwrap();
        assert!(
            fit.confidence < 0.9,
            "flat data reported confidence {}",
            fit.confidence
        );
        assert!(
            fit.log2_ptr.abs() < 0.2,
            "flat data implied log2 PTR {}",
            fit.log2_ptr
        );
    }

    #[test]
    fn too_few_windows_returns_none() {
        assert!(grid_search(&[0, 100], &[1.0, 2.0], 1_000, 16, None).is_none());
        assert!(grid_search(&[0, 100, 200], &[1.0, f64::NAN, 2.0], 1_000, 16, None).is_none());
    }

    #[test]
    fn prior_breaks_a_tie_toward_the_annotation() {
        let g = 1_000_000u64;
        // Two shallow candidate origins; the prior should decide.
        let (pos, mut val) = v_profile(g, 250_000, 0.3, 60);
        for (i, v) in val.iter_mut().enumerate() {
            *v += ((i % 5) as f64 - 2.0) * 0.02;
        }
        let free = grid_search(&pos, &val, g, 128, None).unwrap();
        let primed = grid_search(
            &pos,
            &val,
            g,
            128,
            Some(OriPrior {
                position: 250_000,
                sd: 20_000.0,
            }),
        )
        .unwrap();
        let d_free = circular_distance(free.ori, 250_000, g);
        let d_primed = circular_distance(primed.ori, 250_000, g);
        assert!(
            d_primed <= d_free,
            "prior moved the fit away: {d_primed} vs {d_free}"
        );
    }

    #[test]
    fn annotations_out_of_range_are_ignored() {
        let mut p = std::env::temp_dir();
        p.push(format!("sk2bgrow-{}-ori.tsv", std::process::id()));
        std::fs::write(&p, "# genome\tpos\n g0 \t100\t0.9\tDoriC\nmissing\t50\n").unwrap();
        let ann = read_annotations(&p).unwrap();
        assert_eq!(ann.len(), 2);
        assert_eq!(ann[0].genome, "g0");
        assert_eq!(ann[0].source, "DoriC");
        assert_eq!(
            ann[1].confidence, 1.0,
            "missing confidence should default to 1.0"
        );
        std::fs::remove_file(p).ok();
    }
}
