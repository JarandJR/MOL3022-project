use clap::{Arg, Command};
use de_novo_motif_discover::{
    methods::SupportedMethods,
    motif::{em, gibbs_sampling::GibbsSampling, seq_logo::generate_sequence_logo},
    parser::parse,
};

fn main() {
    let matches = parse_args();
    let method =
        Into::<SupportedMethods>::into(matches.get_one::<String>("method").unwrap().as_str());
    let kmer = *matches.get_one::<usize>("kmer_length").unwrap();
    let path = matches.get_one::<String>("path").unwrap();
    let seqlogo = matches.get_flag("seqlogo");
    let debug = matches.get_flag("debug");
    let treshold = *matches.get_one::<f64>("treshold").unwrap();
    let max_iter = *matches.get_one::<usize>("max_iter").unwrap();

    println!("File path: {}", path);
    println!("Method: {:?}", method);
    println!("Motif length: {}", kmer);
    println!("Treshold for convergence: {}", treshold);
    println!("Maximum number of iterations: {}", max_iter);
    println!("Generate sequence logo: {}", seqlogo);
    println!("Debug: {}\n", debug);

    let (seqs, nc) = parse(path);
    let start = std::time::Instant::now();
    let pwm = match &method {
        &SupportedMethods::Gibbs => {
            GibbsSampling::new(nc, kmer, &seqs).discover_motif(max_iter, treshold, debug)
        }
        &SupportedMethods::EM => em::discover_motif(&seqs, kmer, max_iter, treshold, debug),
        _ => panic!("Unsupported method"),
    };
    let duration = start.elapsed();
    println!("Time taken: {:?}", duration);

    if seqlogo {
        println!("{:?}", pwm);
        generate_sequence_logo(&pwm, &method);
    }
}

fn parse_args() -> clap::ArgMatches {
    let matches = Command::new("Motif Discovery")
        .version("1.0")
        .author("Jarand Romestrand")
        .about("Implements Gibbs Sampler and Expectation Maximization for motif discovery")
        .arg(
            Arg::new("method")
                .short('m')
                .long("method")
                .required(true)
                .value_parser(["gibbs", "em"])
                .help("Choose the algorithm: Gibbs sampler or expectation-maximization"),
        )
        .arg(
            Arg::new("path")
                .short('p')
                .long("path")
                .required(true)
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
        .arg(
            Arg::new("treshold")
                .long("treshold")
                .short('t')
                .default_value("0.005")
                .value_parser(clap::value_parser!(f64))
                .help("Treshold for convergence (default: 0.005)"),
        )
        .arg(
            Arg::new("max_iter")
                .long("max_iter")
                .short('i')
                .default_value("100")
                .value_parser(clap::value_parser!(usize))
                .help("Max iterations for motif discover (default: 100)"),
        )
        .arg(
            Arg::new("debug")
                .long("debug")
                .short('d')
                .action(clap::ArgAction::SetTrue)
                .help("Debug when running"),
        )
        .get_matches();
    matches
}
