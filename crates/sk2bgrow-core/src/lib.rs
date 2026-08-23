//! # sk2bgrow-core
//!
//! Deterministic 2bRAD/TGT anchors as a sketch for bacterial PTR estimation.
//!
//! The layering follows the architecture doc: this crate owns everything from
//! reference sequence to a per-anchor count table. Statistical modelling —
//! ZTP/NB window rates, GC correction, the V-shape fit, cross-enzyme fusion —
//! lives in `python/sk2bgrow/`, and the two halves meet at a file rather than an
//! FFI boundary so they can iterate independently.
//!
//! ```text
//!   genomes/*.fna ──digest──> TGT v2 ──anchor_db──> db/          (offline, once)
//!                                                    │
//!   reads *.fq.gz ────────count───────────────────> counts.tsv   (online, per sample)
//!                                                    │
//!                          python/sk2bgrow ─────────> output.tsv (PTR + QC)
//! ```
//!
//! ## Module map
//!
//! | module | role | reuse source |
//! |---|---|---|
//! | [`seq`] | IUPAC / 2-bit / Hamming primitives | `bsyn::seq` |
//! | [`enzyme`] | the 16 Type IIB panel | `bsyn::enzyme` |
//! | [`digest`] | in silico digestion | `bsyn::digest` |
//! | [`tgt`] | TGT v2 read/write | `bsyn::tgt` |
//! | [`anchor_db`] | anchor library + uniqueness masking | Pilea `sketch.py` (index side) |
//! | [`window`] | equal-anchor adaptive windows | Pilea 25 kb windows |
//! | [`count`] | reads -> anchors, <=2 mismatch | Syn2bANI `tag_matcher.rs` |
//! | [`em`] | shared-anchor reassignment | Pilea/sylph containment EM |
//! | [`ori`] | ori annotation + coarse grid search | new |
//! | [`scaffold`] | MAG contig ordering | `bsyn scaffold` |

pub mod anchor_db;
pub mod count;
pub mod digest;
pub mod em;
pub mod enzyme;
pub mod error;
pub mod fasta;
pub mod ori;
pub mod scaffold;
pub mod seq;
pub mod tgt;
pub mod window;

pub use error::{Result, Sk2bError};

/// Crate version, written into every anchor database manifest.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Convenience re-exports for the CLI and for downstream crates.
pub mod prelude {
    pub use crate::anchor_db::{Anchor, AnchorDb, BuildParams, GenomeMeta};
    pub use crate::count::{AnchorIndex, CountMode, CountStats, MatchConfig};
    pub use crate::digest::DigestConfig;
    pub use crate::em::{EmConfig, EmResult};
    pub use crate::enzyme::{Enzyme, EnzymeSet, Motif, Pattern, PANEL};
    pub use crate::error::{Result, Sk2bError};
    pub use crate::ori::{OriFit, OriPrior};
    pub use crate::tgt::{ContigKind, ContigMeta, Tgt, TgtRecord};
    pub use crate::window::{Window, WindowPolicy};
}
