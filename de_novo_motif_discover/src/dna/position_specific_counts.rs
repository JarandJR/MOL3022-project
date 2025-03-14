use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};

use super::{alignment_matrix::AlignmentMatrix, nucleotide::Nucleotide};

pub struct PositionFrequencyMatrix {
    psc: Vec<Vec<usize>>,
}

impl PositionFrequencyMatrix {
    pub fn new(kmer: usize, alnmtx: &AlignmentMatrix) -> Self {
        let psc =
            alnmtx
                .iter()
                .enumerate()
                .fold(vec![vec![0; kmer]; 4], |mut acc, (idx_seq, i)| {
                    alnmtx
                        .motif(idx_seq)
                        .chars()
                        .enumerate()
                        .for_each(|(i, c)| {
                            if c != 'N' {
                                acc[Into::<usize>::into(Nucleotide::from(c))][i] += 1
                            }
                        });
                    acc
                });
        Self { psc }
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Vec<usize>> {
        self.psc.iter()
    }

    pub fn sum_matrix(&self) -> Vec<usize> {
        self.psc.iter().map(|counts| counts.iter().sum()).collect()
    }

    pub fn add_motif(&mut self, m: &str, x: isize) {
        m.chars().enumerate().for_each(|(i, c)| {
            if c != 'N' {
                self[c][i] = (self[c][i] as isize + x) as usize
            }
        });
    }
}

impl Index<Nucleotide> for PositionFrequencyMatrix {
    type Output = Vec<usize>;

    fn index(&self, index: Nucleotide) -> &Self::Output {
        &self.psc[Into::<usize>::into(index)]
    }
}

impl IndexMut<Nucleotide> for PositionFrequencyMatrix {
    fn index_mut(&mut self, index: Nucleotide) -> &mut Self::Output {
        &mut self.psc[Into::<usize>::into(index)]
    }
}

impl Index<char> for PositionFrequencyMatrix {
    type Output = Vec<usize>;

    fn index(&self, index: char) -> &Self::Output {
        &self[Into::<Nucleotide>::into(index)]
    }
}

impl IndexMut<char> for PositionFrequencyMatrix {
    fn index_mut(&mut self, index: char) -> &mut Self::Output {
        &mut self[Into::<Nucleotide>::into(index)]
    }
}

impl Index<usize> for PositionFrequencyMatrix {
    type Output = Vec<usize>;

    fn index(&self, index: usize) -> &Self::Output {
        &&self[Into::<Nucleotide>::into(index)]
    }
}

impl IndexMut<usize> for PositionFrequencyMatrix {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self[Into::<Nucleotide>::into(index)]
    }
}

impl Debug for PositionFrequencyMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        self.psc
            .iter()
            .enumerate()
            .map(|(nuc, counts)| write!(f, "{}: {:?}\n", Nucleotide::from(nuc), counts))
            .collect()
    }
}
