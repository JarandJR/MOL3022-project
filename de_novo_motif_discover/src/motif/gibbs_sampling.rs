use crate::{
    dna::{alignment_matrix::AlignmentMatrix, nucleotide::NucletideCounts},
    parser::Sequence,
};

use super::pwm::PositionWeightMatrix;

type PositionSpecificCounts = [Vec<usize>; 4];

pub struct GibbsSampling<'a> {
    // Position-Specific Count Matrix
    pscm: PositionSpecificCounts,
    // Alignment matrix
    alnmtx: AlignmentMatrix<'a>,
    // Background count
    bc: NucletideCounts,
    // Nucleotide count
    nc: NucletideCounts,
    // Sequences
    seqs: &'a Vec<Sequence>,
}

impl<'a> GibbsSampling<'a> {
    pub fn new(nc: NucletideCounts, k: usize, seqs: &'a Vec<Sequence>) -> Self {
        const REPEAT_VALUE: Vec<usize> = Vec::new();
        let pscm = [REPEAT_VALUE; 4];
        let alnmtx = AlignmentMatrix::new_random(k, seqs);
        let bc = NucletideCounts::default();
        Self {
            pscm,
            alnmtx,
            bc,
            nc,
            seqs,
        }
    }

    pub fn discover(&mut self) -> &Self {
        //seqs.iter().for_each(|s| println!("{}", s));
        let max_iterations = 1; //_000_000;
        for it in 0..max_iterations {
            println!("Iteration: {it}");
            dbg!(&self.alnmtx);
        }
        self
    }

    pub fn pwm(&self) -> PositionWeightMatrix {
        // Testdata
        PositionWeightMatrix::default()
    }
}
