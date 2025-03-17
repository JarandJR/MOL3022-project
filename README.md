# MOL3022-Project
Bioinformatics - Method Oriented Project

## Overview
This project focuses on predicting **de novo transcription factor binding sites**. It aims to provide a tool that can find potential motifs in DNA sequence data that are likely to represent transcription factor binding sites, which is important in regulating gene expression.

### **Purpose**

The purpose of this project is to implement two algorithms for **de novo motif discovery**:
- **Gibbs Sampler** (stochastic algorithm)
- **Expectation-Maximization (EM)** (deterministic algorithm)

These algorithms are used to identify recurring patterns (motifs) in a set of DNA sequences, which could represent the binding sites of transcription factors.

## Availability
This program is available on all major platforms that support Rust and Python, including:
- Windows
- macOS
- Linux


## Requirements
To run the program, some requirements are needed.
### System requirements
- Rust
- Python

Ensure that you have both Rust and Python installed and properly configured in your environment to run the program successfully.

### Pyhton dependecies for sequence logo
The following Pyhton libraries are required for generation of sequence logo of the discovered motif.
- numpy
- pandas
- logomaker
- matplotlib

You can use the command below to install these requirements:
```bash
pip install numpy pandas logomaker matplotlib
```
### Rust dependecies
The rust dependecies will automatically be downloaded when running the program. These dependecies are as follows:
- clap (for arguments parsing)
- rand (for random number generation)
- serde_json (for passing the pwm to pyhton for sequence logo generation)

## Running the program
To execute the Motif Discovery tool, follow these steps:

### 1. Clone the Repository and Navigate into the Project Directory
First, clone the repository to your local machine:
```bash
git clone https://github.com/JarandJR/MOL3022-project.git
```
Then, change into the project directory:
```bash
cd de_novo_motif_discover
```
### 2. Build and Run the Program
Once inside the project directory, you can run the program using the following command format:
```bash
cargo run --release -- [OPTIONS]
```

### Command-Line Arguments

| Option            | Short | Required | Default  | Description |
|------------------|------|----------|---------|-------------|
| `--method`       | `-m` | ✅ Yes   | N/A     | Choose the algorithm: `gibbs` (Gibbs Sampler) or `em` (Expectation-Maximization). |
| `--path`         | `-p` | ✅ Yes   | N/A     | Path to the sequence data file. |
| `--kmer_length`  | `-k` | ❌ No    | `8`     | Length of the motif to discover. |
| `--logo`         | `-l` | ❌ No    | `false` | Generate a sequence logo. |
| `--treshold`     | `-t` | ❌ No    | `0.005` | Threshold for convergence. |
| `--max_iter`     | `-i` | ❌ No    | `100`   | Maximum iterations for motif discovery. |
| `--debug`        | `-d` | ❌ No    | `false` | Enable debug mode. |
--------------------------------------
#### Other Options
- `--help` or `-h`: Display help information and uage details of the program.  

## Usage
The Motif Discovery tool supports two main algorithms: Gibbs Sampler and Expectation-Maximization (EM). The following examples show how to use the tool with different options on the provided dataset in the `data` folder.

### 1. Running with Gibbs Sampler
The Gibbs Sampler algorithm is used to discover motifs in sequence data stochastically. Here's an example of how to run the program with this method:
```bash
cargo run --release -- -m gibbs -p ../data/CREB1_K562_ChIPseq_79_544.fasta -k 10 -l -t 0.01 -i 200 -d
```
Breakdown:
- `-m` gibbs: Use the Gibbs Sampler algorithm.
- `-p` ../data/CREB1_K562_ChIPseq_79_544.fasta: Path to the sequence data file.
- `-k` 10: Length of the motif to discover (10).
- `-l`: Generate a sequence logo.
- `-t` 0.001: Threshold for convergence (0.001).
- `-i` 200: Maximum iterations (200).
- `-d`: Enable debug mode.

### 2. Running with Expectation-Maximization (EM)
The EM algorithm is used to discover motifs in sequence data deterministic. Here's an example of how to run the program with this method:

```bash
cargo run --release -- -m em -p ../data/CREB1_K562_ChIPseq_79_544.fasta
```
Breakdown:
- `-m` em: Use the Expectation-Maximization algorithm.
- `-p` ../data/CREB1_K562_ChIPseq_79_544.fasta: Path to the sequence data file.
- `-l`: Generate a sequence logo.
- `-i` 200: Maximum iterations (200).
- `-d`: Enable debug mode.

### 3. Running with Default Parameters
If you prefer to run the program with the default parameters, you can omit optional arguments:
```bash
cargo run --release -- -m [METHOD] -p [PATH TO SEQUENCE DATA]
```

In this case, it will:
- Use the default motif length (8).
- Use the default threshold (0.005).
- Use the default maximum iterations (100).

### 4. Displaying Help Information
To get a quick summary of available options and arguments, run:

```bash
cargo run --release -- -h
```

This will show the help information, including all options and their descriptions.

## **Test Dataset**

The test dataset used for development, which is available in the `data` folder, was sourced from the [ENCODE Project](https://www.encodeproject.org/) under **ChIP-seq experiments**.

### **Dataset Details:**
- **Species:** Homo sapiens (Human)
- **Target:** CREB1 (Transcription Factor)
- **Sample:** K562 (A human leukemia cell line, derived from chronic myelogenous leukemia (CML) cells)

This dataset contains sequencing data that can be used for testing the motif discovery tool and evaluating the algorithm's performance. It is provided in the `data` folder, and you are welcome to use it for your own testing and experimentation.

You can access the dataset for further details directly on the [ENCODE Project page](https://www.encodeproject.org/experiments/ENCSR000BSO/).

## Authors
This project was created by:
- Jarand Romestrand – https://github.com/JarandJR – Developer of the Gibbs Sampler implementation
- Vetle Nordang - https://github.com/VetleNordang - Developer of the EM implementation
