from Bio import SeqIO
import random
from tqdm import tqdm

def read_fasta_file(file_path):
    return [str(record.seq) for record in SeqIO.parse(file_path, "fasta")]

def initialize_pwm(sequences, motif_len):

    bases = ["A", "C", "G", "T"]
    pwm = {b: [0.25] * motif_len for b in bases}
    return pwm

def e_step(sequences, pwm, background, motif_len):
    Z = []
    for seq in sequences:
        scores = []
        seq_len = len(seq)
        
        for start in range(seq_len - motif_len + 1):
            window = seq[start:start + motif_len]
            
            if 'N' in window:
                scores.append(0.0)
                continue
                
            # Calculate motif probability
            motif_prob = 1.0
            for j in range(motif_len):
                base = window[j] # base means the nucleotide at position j in the window
                motif_prob *= pwm[base][j]
            
            # Calculate background probability for non-motif positions
            bg_prob = 1.0
            for pos in range(seq_len):
                if pos < start or pos >= start + motif_len:
                    base = seq[pos]
                    if base in background:  # Only multiply if base is in background
                        bg_prob *= background[base]
                    # Skip 'N' characters in background calculation
            
            # Combine both probabilities
            scores.append(motif_prob * bg_prob)
            
        total = sum(scores)
        if total == 0:
            Z.append([1.0 / len(scores)] * len(scores))
        else:
            Z.append([score / total for score in scores])
    
    return Z

def m_step(sequences, Z, motif_len, pseudocount=1e-3):
    pwm = {b: [pseudocount] * motif_len for b in 'ACGT'}

    for seq, probs in zip(sequences, Z):
        for start, weight in enumerate(probs):
            window = seq[start:start + motif_len]
            if 'N' in window:
                continue

            for j in range(motif_len):
                base = window[j]
                pwm[base][j] += weight
    # Normalize each column
    for j in range(motif_len):
        col_sum = sum(pwm[b][j] for b in 'ACGT')
        for b in 'ACGT':
            pwm[b][j] /= col_sum
    return pwm

def run_em(sequences, motif_len, max_iters=100, tol=1e-4):
    pwm = initialize_pwm(sequences, motif_len)

    pbar = tqdm(total=max_iters, desc="EM Algorithm")
    background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    iterations_completed = 0
    for iteration in range(max_iters):
        Z = e_step(sequences, pwm, background, motif_len)
        new_pwm = m_step(sequences, Z, motif_len)
        delta = max(abs(new_pwm[b][j] - pwm[b][j]) for b in 'ACGT' for j in range(motif_len))
        pwm = new_pwm

        iterations_completed += 1
        pbar.update(1)
        pbar.set_postfix({"delta": f"{delta:.6f}", "converged": delta < tol})

        if delta < tol:
            break

    pbar.set_postfix({"delta": f"{delta:.6f}", "converged": delta < tol})
    print(f"EM completed after {iterations_completed} iterations with final delta={delta:.6f}")
    return pwm, Z

def extract_motif_sites(sequences, Z, motif_len):
    sites = []
    for seq, probs in zip(sequences, Z):
        start = probs.index(max(probs))
        sites.append(seq[start:start+motif_len])
    return sites


if __name__ == "__main__":
    file_path = "../../data/CREB1_K562_ChIPseq_79_544.fasta"
    print("Reading FASTA file...")
    sequences = read_fasta_file(file_path)
    print(f"Loaded {len(sequences)} sequences")

    motif_len = 8
    print(f"Starting EM algorithm to find motifs of length {motif_len}...")


    pwm, Z = run_em(sequences, motif_len)
    sites = extract_motif_sites(sequences, Z, motif_len)

    print(f"\nDiscovered {len(sites)} motif sites")
    print("Sample motifs discovered:")
    for i, site in enumerate(sites[:10]):  # Print first 10 motifs
        print(f"Motif {i+1}: {site}")