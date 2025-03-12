use std::{
    fmt::Debug,
    ops::{Index, IndexMut},
};

pub struct NucletideCounts {
    nc: [usize; 5],
}

impl NucletideCounts {
    fn i(c: char) -> usize {
        match c {
            'A' => 0,
            'T' => 1,
            'G' => 2,
            'C' => 3,
            'N' => 4,
            _ => panic!("Nucleotide error"),
        }
    }

    fn c(i: usize) -> char {
        match i {
            0 => 'A',
            1 => 'T',
            2 => 'G',
            3 => 'C',
            4 => 'N',
            _ => panic!("Nucleotide error"),
        }
    }
}

impl Default for NucletideCounts {
    fn default() -> Self {
        Self { nc: [0; 5] }
    }
}

impl Index<char> for NucletideCounts {
    type Output = usize;

    fn index(&self, index: char) -> &Self::Output {
        &self.nc[Self::i(index)]
    }
}

impl IndexMut<char> for NucletideCounts {
    fn index_mut(&mut self, index: char) -> &mut usize {
        &mut self.nc[Self::i(index)]
    }
}

impl Debug for NucletideCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.nc
            .iter()
            .enumerate()
            .map(|(i, nc)| write!(f, "{}: {nc} ", Self::c(i)))
            .collect()
    }
}
