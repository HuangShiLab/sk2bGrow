//! Error type shared across the core crate.

use std::path::PathBuf;

pub type Result<T> = std::result::Result<T, Sk2bError>;

#[derive(Debug, thiserror::Error)]
pub enum Sk2bError {
    #[error("I/O error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("malformed input {path}: {msg}")]
    Format { path: PathBuf, msg: String },

    #[error("anchor database error: {0}")]
    Db(String),

    #[error(transparent)]
    Enzyme(#[from] crate::enzyme::EnzymeError),

    #[error("serialisation error: {0}")]
    Serde(String),

    #[error("{0}")]
    Config(String),
}

impl From<bincode::Error> for Sk2bError {
    fn from(e: bincode::Error) -> Self {
        Sk2bError::Serde(e.to_string())
    }
}

impl From<serde_json::Error> for Sk2bError {
    fn from(e: serde_json::Error) -> Self {
        Sk2bError::Serde(e.to_string())
    }
}
