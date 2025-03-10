use std::{fs::read_to_string, path::Path};

use fasta::FastaSeq;

use crate::filetypes::SupportedFiletypes;

pub mod fasta;

pub type Sequence = String;

pub fn parse(path: &str) -> Vec<Sequence> {
    let file_type = Into::<SupportedFiletypes>::into(
        Path::new(path)
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("unknown"),
    );
    let contents = read_to_string(path).unwrap();
    match file_type {
        SupportedFiletypes::Fasta => FastaSeq::sequences(contents),
        _ => panic!("Unsupported file type"),
    }
}
