use std::fmt::Debug;

use rand::random_range;

use crate::parser::Sequence;

pub struct AlignmentMatrix<'a> {
    alin: Vec<usize>,
    seqs: &'a Vec<Sequence>,
    kmer: usize,
}

impl<'a> AlignmentMatrix<'a> {
    pub fn new_random(kmer: usize, seqs: &'a Vec<Sequence>) -> Self {
        let alin = Self::initialize(kmer, seqs);
        Self { alin, seqs, kmer }
    }

    fn initialize(kmer: usize, seqs: &'a Vec<Sequence>) -> Vec<usize> {
        seqs.iter()
            .map(|s| {
                let max = s.len() - kmer;
                random_range(0..=max)
            })
            .collect()
    }
}

// Todo implement actual alignement debug
impl<'a> Debug for AlignmentMatrix<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        self.alin
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let max = self.seqs[i].len() - self.kmer;
                write!(
                    f,
                    "{}{}{}{}\n",
                    " ".repeat(max - a),
                    &self.seqs[i][0..*a].to_lowercase(),
                    &self.seqs[i][*a..a + self.kmer],
                    &self.seqs[i][a + self.kmer..].to_lowercase(),
                )
            })
            .collect()
    }
}
