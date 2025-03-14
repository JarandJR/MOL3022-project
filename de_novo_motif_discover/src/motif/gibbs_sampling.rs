use rand::{rngs::ThreadRng, Rng};

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
    // Background Frequency
    bgf: NucletideCounts,
    // Sequences
    seqs: &'a Vec<Sequence>,
}

impl<'a> GibbsSampling<'a> {
    pub fn new(nc: NucletideCounts, k: usize, seqs: &'a Vec<Sequence>) -> Self {
        let alnmtx = AlignmentMatrix::new_random(k, seqs);
        let pfm = PositionFrequencyMatrix::new(k, &alnmtx);
        let bgf = pfm
            .sum_matrix()
            .iter()
            .zip(nc.iter())
            .map(|(nc, tot)| tot - nc)
            .collect::<Vec<_>>()
            .into();
        Self {
            pfm,
            alnmtx,
            bgf,
            seqs,
        }
    }

    pub fn discover(
        &mut self,
        max_iterations: usize,
        treshold: f64,
        debug: bool,
        print_align: bool,
    ) -> PositionWeightMatrix {
        let mut max_score = self.score();
        let mut prev_scores = Vec::new();
        let mut max_pwm = self.pwm().collect::<Vec<Vec<_>>>().into();
        let mut rng = rand::rng();

        debug.then(|| println!("Start score: {}", max_score));

        for it in 0..max_iterations {
            debug.then(|| println!("\nIteration: {}", it));
            let mut order = (0..self.seqs.len()).collect::<Vec<_>>();
            while !order.is_empty() {
                let idx_order = rng.random_range(0..order.len());
                let idx_seq = order.remove(idx_order);
                let seq = self.seqs[idx_seq].as_str();

                // Take out the current motif for prediction
                let motif = self.alnmtx.motif(idx_seq);
                self.bgf.add_motif(motif, 1);
                self.pfm.add_motif(motif, -1);

                let ppm = self.position_probability_matrix().collect::<Vec<_>>();
                let bgp = self.background_probability().collect::<Vec<_>>();

                // Calculate the probability-score for all possible motif starting positions
                let ps = (0..=seq.len() - self.kmer())
                    // Motifs
                    .map(|start_pos| &seq[start_pos..start_pos + self.kmer()])
                    // Probability score
                    .map(|m| {
                        m.chars().enumerate().fold(1., |acc, (i, c)| {
                            if c != 'N' {
                                let probability = ppm[Into::<usize>::into(Nucleotide::from(c))][i]
                                    / bgp[Into::<usize>::into(Nucleotide::from(c))];
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
                let new_start_pos = self.pick_random_from_normalized(&mut rng, norm_ps);

                // Update alignment matrix and add new motif to frequency matrix
                self.alnmtx.update_pos(idx_seq, new_start_pos);
                let new_motif = self.alnmtx.motif(idx_seq);
                self.bgf.add_motif(new_motif, -1);
                self.pfm.add_motif(new_motif, 1);
            }
            let new_score = self.score();
            debug.then(|| println!("Score: {}", new_score));
            if max_score < new_score {
                max_score = new_score;
                max_pwm = self.pwm().collect::<Vec<Vec<_>>>().into();
                continue;
            }
            prev_scores.push(new_score);
            let mean = prev_scores.iter().sum::<f64>() / prev_scores.len() as f64;
            let relative_imp = (max_score - mean).abs() / mean;
            if relative_imp < treshold {
                debug.then(|| println!("Converged"));
                break;
            }
        }
        print_align.then(|| println!("{:?}\n", self.alnmtx));
        debug.then(|| println!("Final score: {}", max_score));
        max_pwm
    }

    fn pick_random_from_normalized(&self, rng: &mut ThreadRng, norm: Vec<f64>) -> usize {
        let mut cumulative_sum = 0.0;
        let roll = rng.random_range(0.0..1.);
        norm.iter()
            .enumerate()
            .find_map(|(i, ps)| {
                cumulative_sum += ps;
                cumulative_sum.ge(&roll).then_some(i)
            })
            .unwrap_or(norm.len() - 1)
    }

    fn score(&self) -> f64 {
        self.pfm
            .iter()
            .zip(self.pwm())
            .map(|(freqs, weights)| {
                freqs
                    .iter()
                    .zip(weights)
                    .map(|(c, w)| *c as f64 * w)
                    .sum::<f64>()
            })
            .sum()
    }

    fn background_probability(&self) -> impl Iterator<Item = f64> + '_ {
        let tot = self.bgf.iter().sum::<usize>() as f64;
        self.bgf.iter().map(move |c| *c as f64 / tot)
    }

    fn kmer(&self) -> usize {
        self.alnmtx.kmer
    }

    fn position_probability_matrix(&self) -> impl Iterator<Item = Vec<f64>> + '_ {
        let pseudocount = 0.0001;
        let total = self.seqs.len() as f64 + 4. * pseudocount;
        self.pfm.iter().map(move |nuc_freqs| {
            nuc_freqs
                .iter()
                .map(|freq| (*freq as f64 + pseudocount) / total)
                .collect::<Vec<_>>()
        })
    }

    fn pwm(&self) -> impl Iterator<Item = Vec<f64>> + '_ {
        self.position_probability_matrix()
            .zip(self.background_probability())
            .map(|(probabilities, bgp)| {
                probabilities
                    .iter()
                    .map(|p| {
                        let pos_weight = p / bgp;
                        pos_weight.log2()
                    })
                    .collect::<Vec<_>>()
            })
    }
}
