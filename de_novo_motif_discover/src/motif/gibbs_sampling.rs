use rand::Rng;

use crate::{
    dna::{
        alignment_matrix::AlignmentMatrix, nucleotide::Nucleotide,
        nucleotide_count::NucletideCounts, position_specific_counts::PositionFrequencyMatrix,
    },
    parser::Sequence,
};

use super::pwm::PositionWeightMatrix;

pub struct GibbsSampling<'a> {
    // Position Frequency Matrix
    pfm: PositionFrequencyMatrix,
    // Alignment matrix
    alnmtx: AlignmentMatrix<'a>,
    // Background count
    bc: NucletideCounts,
    // Sequences
    seqs: &'a Vec<Sequence>,
}

impl<'a> GibbsSampling<'a> {
    pub fn new(nc: NucletideCounts, k: usize, seqs: &'a Vec<Sequence>) -> Self {
        let alnmtx = AlignmentMatrix::new_random(k, seqs);
        let psc = PositionFrequencyMatrix::new(k, &alnmtx);
        let bc = psc
            .sum_matrix()
            .iter()
            .zip(nc.iter())
            .map(|(nc, tot)| tot - nc)
            .collect::<Vec<_>>()
            .into();
        Self {
            pfm: psc,
            alnmtx,
            bc,
            seqs,
        }
    }

    pub fn discover(&mut self, max_iterations: usize) -> &Self {
        let mut prev_score = self.score();
        let mut rng = rand::rng();
        for it in 0..max_iterations {
            println!("iteration: {}", it);
            let mut order = (0..self.seqs.len()).collect::<Vec<_>>();
            while !order.is_empty() {
                let idx_order = rng.random_range(0..order.len());
                let idx_seq = order.remove(idx_order);
                let seq = &self.seqs[idx_seq];

                // Take out the current motif for prediction
                let motif = self.alnmtx.motif(idx_seq);
                self.bc.add_motif(motif, 1);
                self.pfm.add_motif(motif, -1);

                let ppm = self.position_probability_matrix();
                let bgf = self.background_frequency().collect::<Vec<_>>();

                // Calculate the probability-score for all possible motif starting positions
                let ps = (0..=seq.len() - self.kmer())
                    // Motifs
                    .map(|start_pos| &seq[start_pos..start_pos + self.kmer()])
                    // Probability score
                    .map(|m| {
                        m.chars().enumerate().fold(1., |acc, (i, c)| {
                            if c != 'N' {
                                let probability = ppm[Into::<usize>::into(Nucleotide::from(c))][i]
                                    / bgf[Into::<usize>::into(Nucleotide::from(c))];
                                return acc * probability;
                            }
                            acc
                        })
                    })
                    .collect::<Vec<_>>();

                // Normalize the probability-scores
                let sum_ps = ps.iter().sum::<f64>();
                let norm_ps = ps.iter().map(|s| s / sum_ps).collect::<Vec<_>>();

                // Stochastically determine a new starting position for the motif
                let mut cumulative_sum = 0.0;
                let roll = rng.random_range(0.0..1.);
                let new_start_pos = norm_ps
                    .iter()
                    .enumerate()
                    .find_map(|(i, ps)| {
                        cumulative_sum += ps;
                        cumulative_sum.ge(&roll).then_some(i)
                    })
                    .unwrap_or(norm_ps.len() - 1);

                // Update alignment matrix and add new motif to frequency matrix
                self.alnmtx.update_pos(idx_seq, new_start_pos);
                let new_motif = self.alnmtx.motif(idx_seq);
                self.bc.add_motif(new_motif, -1);
                self.pfm.add_motif(new_motif, 1);
            }
            // TODO, some treshold thing, huge boost in comp time
            let new_score = dbg!(self.score());
            if new_score < prev_score {
                //break;
            }
            prev_score = new_score;
        }
        self
    }

    fn score(&self) -> f64 {
        let pwm = self.pwm();
        self.pfm
            .iter()
            .zip(pwm.iter())
            .map(|(counts, weights)| {
                counts
                    .iter()
                    .zip(weights)
                    .map(|(c, w)| *c as f64 * w)
                    .sum::<f64>()
            })
            .sum()
    }

    fn background_frequency(&self) -> impl Iterator<Item = f64> + '_ {
        let tot = self.bc.iter().sum::<usize>() as f64;
        self.bc.iter().map(move |c| *c as f64 / tot)
    }

    fn kmer(&self) -> usize {
        self.alnmtx.kmer
    }

    fn position_probability_matrix(&self) -> Vec<Vec<f64>> {
        let pseudocount = 0.0001;
        let total = self.seqs.len() as f64 + 4. * pseudocount;
        self.pfm
            .iter()
            .map(|nuc_counts| {
                nuc_counts
                    .iter()
                    .map(|count| (*count as f64 + pseudocount) / total)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<Vec<_>>>()
    }

    pub fn pwm(&self) -> PositionWeightMatrix {
        let pseudocount = 0.0001;
        let total = self.seqs.len() as f64 + 4. * pseudocount;
        self.pfm
            .iter()
            .zip(self.background_frequency())
            .map(|(nuc_counts, bgf)| {
                nuc_counts
                    .iter()
                    .map(|count| {
                        let pos_probability = (*count as f64 + pseudocount) / total;
                        let pos_weight = pos_probability / bgf;
                        pos_weight.log2()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<Vec<_>>>()
            .into()
    }
}
