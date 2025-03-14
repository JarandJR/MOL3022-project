use std::fmt::Debug;

use rand::random_range;

use crate::parser::Sequence;

pub struct AlignmentMatrix<'a> {
    alin: Vec<usize>,
    seqs: &'a Vec<Sequence>,
    pub kmer: usize,
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

    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.alin.iter()
    }

    pub fn motif(&self, idx_seq: usize) -> &str {
        let idx = self.alin[idx_seq];
        &self.seqs[idx_seq][idx..idx + self.kmer]
    }

    pub fn update_pos(&mut self, idx_seq: usize, new_pos: usize) {
        self.alin[idx_seq] = new_pos;
    }
}

impl<'a> Debug for AlignmentMatrix<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        self.alin
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let max = self.seqs[i].len() - self.kmer;
                let prefix = &self.seqs[i][0..*a].to_lowercase();
                let motif = self.motif(i);
                let suffix = &self.seqs[i][a + self.kmer..].to_lowercase();
                // ANSI escape codes for Blue and bold
                let colored_motif = format!("\x1b[1;34m{}\x1b[0m", motif);
                write!(
                    f,
                    "{}{}{}{}\n",
                    " ".repeat(max - a),
                    prefix,
                    colored_motif,
                    suffix
                )
            })
            .collect()
    }
}
