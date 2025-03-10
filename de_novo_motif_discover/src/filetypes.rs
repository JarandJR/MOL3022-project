#[derive(Debug)]
pub enum SupportedFiletypes {
    Fasta,
    Unknown,
}

impl From<&str> for SupportedFiletypes {
    fn from(value: &str) -> Self {
        match value {
            "fasta" => Self::Fasta,
            _ => Self::Unknown,
        }
    }
}
