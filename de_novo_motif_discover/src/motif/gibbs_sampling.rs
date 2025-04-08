use indicatif::{ProgressBar, ProgressStyle};
use rand::{rngs::ThreadRng, Rng};

use crate::{
    dna::{
        alignment_matrix::AlignmentMatrix, nucleotide::Nucleotide,
        nucleotide_count::NucletideCounts, position_specific_counts::PositionFrequencyMatrix,
    },
    motif::em::generate_consensus_from_pwm,
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
    /// Initializes the frequencies for sampling
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

    pub fn discover_motif(
        &mut self,
        max_iterations: usize,
        convergence_treshold: f64,
        debug: bool,
    ) -> PositionWeightMatrix {
        let mut rng = rand::rng();
        let mut max_score = self.score();
        let mut prev_rel_improvements = Vec::new();
        let mut max_pwm = self.pwm().collect::<Vec<Vec<_>>>().into();

        let pb = self.get_progress_bar(max_iterations);
        pb.set_message(format!("Gibbs Sampler | score: {:.6}", max_score));

        for _ in 0..max_iterations {
            // Random order of all sequences
            let mut order = (0..self.seqs.len()).collect::<Vec<_>>();
            while !order.is_empty() {
                let idx_order = rng.random_range(0..order.len());
                let idx_seq = order.remove(idx_order);
                let seq = self.seqs[idx_seq].as_str();

                // Take out the current motif from the frequency matrix for prediction of new starting position
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
            pb.inc(1);
            // Update the alignment
            if max_score < new_score {
                max_score = new_score;
                max_pwm = self.pwm().collect::<Vec<Vec<_>>>().into();
                pb.set_message(format!(
                    "Gibbs Sampler | max: {:.2} score: {:.2}, converged: ",
                    max_score, new_score
                ));
                continue;
            }
            // Relative improvement with decay
            let decay = (prev_rel_improvements.len() + 1) as f64;
            prev_rel_improvements.push((max_score - new_score) / (max_score * decay));
            let mean_relative_imp =
                prev_rel_improvements.iter().sum::<f64>() / prev_rel_improvements.len() as f64;
            let mean_relative_imp = mean_relative_imp.abs();
            pb.set_message(format!(
                "Gibbs Sampler | max: {:.2} score: {:.2}, converged: {:.4}",
                max_score, new_score, mean_relative_imp
            ));
            // Check for convergence
            if mean_relative_imp < convergence_treshold {
                pb.inc(1);
                pb.set_message(format!(
                    "Gibbs Sampler | max: {:.2} score: {:.2}, converged: {:.4}",
                    max_score, new_score, mean_relative_imp
                ));
                break;
            }
        }
        pb.finish();
        debug.then(|| {
            println!("\n\n===== Motif =====");
            println!("Final score: {}", max_score);
            println!("Position Weight Matrix:");
            println!("{:?}", max_pwm);
            let consensus = generate_consensus_from_pwm(&max_pwm);
            println!("Consensus sequence: {}", consensus);
        });
        max_pwm
    }

    fn get_progress_bar(&self, max_iters: usize) -> ProgressBar {
        let pb = ProgressBar::new(max_iters as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{msg} [{bar:40}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("=> "),
        );
        pb.set_message("EM Algorithm");
        pb
    }

    /// Picks a random variable from a normalized list of data using cumulative sum
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

    /// Scores the current alignment
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

    /// Calculates the background probability for each nucleotide
    fn background_probability(&self) -> impl Iterator<Item = f64> + '_ {
        let tot = self.bgf.iter().sum::<usize>() as f64;
        self.bgf.iter().map(move |c| *c as f64 / tot)
    }

    /// Returns the kmer-length being search for
    fn kmer(&self) -> usize {
        self.alnmtx.kmer
    }

    /// Calculates the Position Probability Matrix for the current alignment
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

    /// Calculates the Position Weight Matrix for the current alignment
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
