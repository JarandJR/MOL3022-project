use hashbrown::HashMap;
use indicatif::{ProgressBar, ProgressStyle};

use crate::dna::nucleotide::Nucleotide;
use crate::parser::Sequence;

use super::pwm::PositionWeightMatrix;

type PWM = HashMap<char, Vec<f64>>;

pub fn discover_motif(
    seqs: &Vec<Sequence>,
    kmer: usize,
    max_iter: usize,
    treshold: f64,
    debug: bool,
) -> PositionWeightMatrix {
    // Run EM on current sequences
    let (pwm, z) = run_em(seqs, kmer, max_iter, treshold);
    // Extract motif sites for this run
    let sites = extract_motif_sites(seqs, &z, kmer);
    let pwm = ['A', 'C', 'G', 'T']
        .into_iter()
        .map(|b| pwm.get(&b).unwrap().to_owned())
        .collect::<Vec<Vec<f64>>>();
    let pwm = pwm.into();
    debug.then(|| {
        println!("\n===== Motif =====");
        println!("Position Weight Matrix:");
        println!("{:?}", pwm);

        // Print sample sites
        println!("\nDiscovered {} sites for this motif", sites.len());
        println!("Sample sites:");
        for (j, site) in sites.iter().take(5).enumerate() {
            println!("Site {}: {}", j + 1, site);
        }
        let consensus = generate_consensus_from_pwm(&pwm);
        println!("Consensus sequence: {}", consensus);
    });
    pwm
}

fn run_em(
    sequences: &[String],
    motif_len: usize,
    max_iters: usize,
    tol: f64,
) -> (PWM, Vec<Vec<f64>>) {
    // Fix: Return only PWM and Z
    // Initialize
    let mut pwm = initialize_pwm(motif_len);

    // Initial background
    let mut background = HashMap::new();
    background.insert('A', 0.25);
    background.insert('C', 0.25);
    background.insert('G', 0.25);
    background.insert('T', 0.25);

    // Setup progress bar
    let pb = ProgressBar::new(max_iters as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{msg} [{bar:40}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=> "),
    );
    pb.set_message("EM Algorithm");

    let mut iterations_completed = 0;
    let mut delta = 0.0;

    for _ in 0..max_iters {
        // Fix: prefix unused variable with _
        // E-step using current PWM and background
        let z = e_step(sequences, &pwm, &background, motif_len);

        // M-step returning updated PWM and background
        let (new_pwm, new_background) = m_step(sequences, &z, motif_len, 1e-3);

        // Calculate delta
        delta = 0.0;
        for &base in &['A', 'C', 'G', 'T'] {
            // Check delta for motif PWM
            for j in 0..motif_len {
                let diff = (new_pwm.get(&base).unwrap()[j] - pwm.get(&base).unwrap()[j]).abs();
                if diff > delta {
                    delta = diff;
                }
            }

            // Check delta for background
            let bg_diff =
                (new_background.get(&base).unwrap() - background.get(&base).unwrap()).abs();
            if bg_diff > delta {
                delta = bg_diff;
            }
        }

        // Update models
        pwm = new_pwm;
        background = new_background;

        // Update progress
        iterations_completed += 1;
        pb.inc(1);
        pb.set_message(format!(
            "EM Algorithm | delta: {:.6}, converged: {}",
            delta,
            delta < tol
        ));

        if delta < tol {
            break;
        }
    }

    pb.finish();
    println!(
        "EM completed after {} iterations with final delta={:.6}",
        iterations_completed, delta
    );
    let z = e_step(sequences, &pwm, &background, motif_len);
    (pwm, z)
}

fn initialize_pwm(motif_len: usize) -> PWM {
    let bases = vec!['A', 'C', 'G', 'T'];
    let mut pwm = HashMap::new();

    for &base in &bases {
        pwm.insert(base, vec![0.25; motif_len]);
    }

    pwm
}

fn e_step(
    sequences: &[String],
    pwm: &PWM,
    background: &HashMap<char, f64>,
    motif_len: usize,
) -> Vec<Vec<f64>> {
    let mut z = Vec::with_capacity(sequences.len());

    for seq in sequences {
        let seq_len = seq.len();
        let mut scores = Vec::new();

        for start in 0..=seq_len.saturating_sub(motif_len) {
            let window = &seq[start..start + motif_len];

            if window.contains('N') {
                scores.push(0.0);
                continue;
            }

            // Calculate motif probability
            let mut motif_prob = 1.0;
            for (j, base) in window.chars().enumerate() {
                motif_prob *= pwm.get(&base).map_or(0.0, |probs| probs[j]);
            }

            // Calculate background probability for non-motif positions
            let mut bg_prob = 1.0;
            for (pos, base) in seq.chars().enumerate() {
                if pos < start || pos >= start + motif_len {
                    if let Some(&prob) = background.get(&base) {
                        bg_prob *= prob;
                    }
                    // Skip 'N' characters in background calculation
                }
            }

            // Combine both probabilities
            scores.push(motif_prob * bg_prob);
        }

        let total: f64 = scores.iter().sum();
        if total == 0.0 {
            let uniform_prob = 1.0 / scores.len() as f64;
            z.push(vec![uniform_prob; scores.len()]);
        } else {
            z.push(scores.iter().map(|&score| score / total).collect());
        }
    }

    z
}

fn m_step(
    sequences: &[String],
    z: &[Vec<f64>],
    motif_len: usize,
    pseudocount: f64,
) -> (PWM, HashMap<char, f64>) {
    // Initialize PWM for the motif
    let mut pwm: PWM = HashMap::new();
    for &base in &['A', 'C', 'G', 'T'] {
        pwm.insert(base, vec![pseudocount; motif_len]);
    }

    // Initialize background counts
    let mut bg_counts: HashMap<char, f64> = HashMap::new();
    for &base in &['A', 'C', 'G', 'T'] {
        bg_counts.insert(base, pseudocount);
    }

    // Count total weighted bases for background calculation
    let mut total_bg_positions: f64 = 0.0;

    // Process each sequence
    for (_, (seq, probs)) in sequences.iter().zip(z.iter()).enumerate() {
        // Fix: prefix unused variable
        // For each possible start position
        for (start, &weight) in probs.iter().enumerate() {
            if start + motif_len > seq.len() {
                continue;
            }

            let window = &seq[start..start + motif_len];
            if window.contains('N') {
                continue;
            }

            // Update motif counts
            for (j, base) in window.chars().enumerate() {
                if let Some(counts) = pwm.get_mut(&base) {
                    counts[j] += weight;
                }
            }

            // Update background counts from non-motif positions
            for (pos, base) in seq.chars().enumerate() {
                // Skip positions inside the motif
                if pos >= start && pos < start + motif_len {
                    continue;
                }

                // Skip 'N' characters
                if base == 'N' {
                    continue;
                }

                // Add to background with weight
                if let Some(count) = bg_counts.get_mut(&base) {
                    *count += weight;
                    total_bg_positions += weight;
                }
            }
        }
    }

    // Normalize motif PWM
    for j in 0..motif_len {
        let col_sum: f64 = ['A', 'C', 'G', 'T']
            .iter()
            .map(|&b| pwm.get(&b).map_or(0.0, |counts| counts[j]))
            .sum();

        for &base in &['A', 'C', 'G', 'T'] {
            if let Some(counts) = pwm.get_mut(&base) {
                counts[j] /= col_sum;
            }
        }
    }

    // Normalize background
    let mut background = HashMap::new();
    for &base in &['A', 'C', 'G', 'T'] {
        background.insert(base, *bg_counts.get(&base).unwrap() / total_bg_positions);
    }

    (pwm, background)
}

fn extract_motif_sites(sequences: &[String], z: &[Vec<f64>], motif_len: usize) -> Vec<String> {
    let mut sites = Vec::new();

    for (seq, probs) in sequences.iter().zip(z.iter()) {
        if let Some(max_idx) = probs
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        {
            let start = max_idx.0;
            if start + motif_len <= seq.len() {
                sites.push(seq[start..start + motif_len].to_string());
            }
        }
    }

    sites
}

fn generate_consensus_from_pwm(pwm: &PositionWeightMatrix) -> String {
    let kmer = pwm[Nucleotide::A].len();
    let mut consensus = String::with_capacity(kmer);

    for pos in 0..kmer {
        let mut max_base = 'N';
        let mut max_prob = 0.0;

        for &base in &['A', 'C', 'G', 'T'] {
            let probs = &pwm[base.into()];
            if probs[pos] > max_prob {
                max_prob = probs[pos];
                max_base = base;
            }
        }

        consensus.push(max_base);
    }

    consensus
}
