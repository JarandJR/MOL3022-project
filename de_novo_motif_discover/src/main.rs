use clap::{Arg, Command};
use de_novo_motif_discover::{
    methods::SupportedMethods,
    motif::{gibbs_sampling::GibbsSampling, seq_logo::generate_sequence_logo},
    parser::parse,
};

fn main() {
    let matches = Command::new("Motif Discovery")
        .version("1.0")
        .author("Your Name")
        .about("Implements Gibbs Sampler and Expectation Maximization for motif discovery")
        .arg(
            Arg::new("method")
                .short('m')
                .long("method")
                //.required(true)
                .default_value("gibbs") // Remove this later
                .value_parser(["gibbs", "em"])
                .help("Choose the algorithm: Gibbs sampler or expectation-maximization"),
        )
        .arg(
            Arg::new("path")
                .short('p')
                .long("path")
                //.required(true)
                .default_value("../data/CREB1_K562_ChIPseq_79_544.fasta") // Remove this later
                .help("Path to the sequence data file"),
        )
        .arg(
            Arg::new("kmer_length")
                .short('k')
                .long("kmer_length")
                .default_value("8")
                .value_parser(clap::value_parser!(usize))
                .help("Length of the motif to discover (default: 8)"),
        )
        .arg(
            Arg::new("seqlogo")
                .long("logo")
                .short('l')
                .action(clap::ArgAction::SetTrue)
                .help("Generate a sequence logo"),
        )
        .get_matches();

    let method =
        Into::<SupportedMethods>::into(matches.get_one::<String>("method").unwrap().as_str());
    let k = *matches.get_one::<usize>("kmer_length").unwrap();
    let path = matches.get_one::<String>("path").unwrap();
    let seqlogo = matches.get_flag("seqlogo");

    println!("Method: {:?}", method);
    println!("Motif length: {}", k);
    println!("File path: {}", path);
    println!("Generate sequence logo: {}", seqlogo);

    let (seqs, nc) = parse(path);
    let start = std::time::Instant::now();
    let pwm = match &method {
        &SupportedMethods::Gibbs => GibbsSampling::new(nc, k, &seqs).discover(100).pwm(),
        &SupportedMethods::EM => panic!("not implemented"),
        _ => panic!("Unsupported method"),
    };
    let duration = start.elapsed();
    println!("Time taken: {:?}", duration);

    if seqlogo {
        println!("{:?}", pwm);
        generate_sequence_logo(&pwm, &method);
    }
}
