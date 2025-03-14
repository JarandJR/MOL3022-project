use super::Sequence;
use crate::dna::nucleotide_count::NucletideCounts;

pub struct FastaSeq;

impl FastaSeq {
    pub fn sequences(contents: String) -> (Vec<Sequence>, NucletideCounts) {
        let mut nc = NucletideCounts::default();
        let seqs = contents
            .lines()
            .skip(1)
            .step_by(2)
            .map(|s| {
                let seq = s.to_owned();
                seq.chars().for_each(|c| nc[c] += 1);
                seq
            })
            .collect();
        (seqs, nc)
    }
}
