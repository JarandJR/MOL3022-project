use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

fn main() -> io::Result<()> {
    let n = 1_000_000;
    let formatted = format_with_underscores(n / 2);
    let input_path = "../CREB1_K562_ChIPseq_31_300_055.fasta";
    let output_path = format!("../CREB1_K562_ChIPseq_{}.fasta", formatted);
    
    let input_file = File::open(input_path)?;
    let reader = BufReader::new(input_file);
    let mut output_file = File::create(output_path.as_str())?;
    
    for (i, line) in reader.lines().enumerate() {
        if i >= n {
            break;
        }
        writeln!(output_file, "{}", line?)?;
    }
    
    println!("Successfully written first {} lines to {}",n, output_path);
    Ok(())
}

fn format_with_underscores(n: usize) -> String {
    let num_str = n.to_string();
    let mut result = String::new();
    let len = num_str.len();
    
    for (i, c) in num_str.chars().enumerate() {
        if i > 0 && (len - i) % 3 == 0 {
            result.push('_');
        }
        result.push(c);
    }
    result
}