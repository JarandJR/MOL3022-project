use std::fmt::Display;

pub enum Nucleotide {
    A,
    C,
    G,
    T,
    N,
}

impl From<char> for Nucleotide {
    fn from(value: char) -> Self {
        match value {
            'A' => Nucleotide::A,
            'C' => Nucleotide::C,
            'G' => Nucleotide::G,
            'T' => Nucleotide::T,
            'N' => Nucleotide::N,
            _ => panic!("Nucleotide error"),
        }
    }
}

impl Into<char> for &Nucleotide {
    fn into(self) -> char {
        match self {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
            Nucleotide::N => 'N',
        }
    }
}

impl From<usize> for Nucleotide {
    fn from(value: usize) -> Self {
        match value {
            0 => Nucleotide::A,
            1 => Nucleotide::C,
            2 => Nucleotide::G,
            3 => Nucleotide::T,
            4 => Nucleotide::N,
            _ => panic!("Nucleotide error"),
        }
    }
}

impl Into<usize> for Nucleotide {
    fn into(self) -> usize {
        match self {
            Nucleotide::A => 0,
            Nucleotide::C => 1,
            Nucleotide::G => 2,
            Nucleotide::T => 3,
            Nucleotide::N => 4,
        }
    }
}

impl Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Into::<char>::into(self))
    }
}
