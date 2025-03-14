use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};

use super::nucleotide::Nucleotide;

pub struct NucletideCounts {
    nc: Vec<usize>,
}

impl NucletideCounts {
    pub fn new(nc: Vec<usize>) -> Self {
        Self { nc }
    }

    pub fn iter(&self) -> std::slice::Iter<'_, usize> {
        self.nc.iter()
    }

    pub fn add_motif(&mut self, m: &str, x: isize) {
        m.chars().for_each(|c| {
            if c != 'N' {
                self[c] = (self[c] as isize + x) as usize
            }
        });
    }
}

impl Default for NucletideCounts {
    fn default() -> Self {
        let nc = vec![0; 5];
        Self::new(nc)
    }
}

impl From<Vec<usize>> for NucletideCounts {
    fn from(value: Vec<usize>) -> Self {
        Self::new(value)
    }
}

impl Index<Nucleotide> for NucletideCounts {
    type Output = usize;

    fn index(&self, index: Nucleotide) -> &Self::Output {
        &self.nc[Into::<usize>::into(index)]
    }
}

impl IndexMut<Nucleotide> for NucletideCounts {
    fn index_mut(&mut self, index: Nucleotide) -> &mut Self::Output {
        &mut self.nc[Into::<usize>::into(index)]
    }
}

impl Index<char> for NucletideCounts {
    type Output = usize;

    fn index(&self, index: char) -> &Self::Output {
        &self[Into::<Nucleotide>::into(index)]
    }
}

impl IndexMut<char> for NucletideCounts {
    fn index_mut(&mut self, index: char) -> &mut usize {
        &mut self[Into::<Nucleotide>::into(index)]
    }
}

impl Index<usize> for NucletideCounts {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &&self[Into::<Nucleotide>::into(index)]
    }
}

impl IndexMut<usize> for NucletideCounts {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self[Into::<Nucleotide>::into(index)]
    }
}

impl Debug for NucletideCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.nc
            .iter()
            .enumerate()
            .map(|(i, nc)| write!(f, "{}: {nc} ", Into::<Nucleotide>::into(i)))
            .collect()
    }
}
