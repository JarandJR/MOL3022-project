use std::fmt;

pub enum SupportedMethods {
    Gibbs,
    EM,
    Unsupported,
}

impl fmt::Debug for SupportedMethods {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Gibbs => write!(f, "Gibbs sampler"),
            Self::EM => write!(f, "Expectation Maximization"),
            Self::Unsupported => write!(f, "Unsupported"),
        }
    }
}

impl From<&str> for SupportedMethods {
    fn from(value: &str) -> Self {
        match value {
            "gibbs" => Self::Gibbs,
            "em" => Self::EM,
            _ => Self::Unsupported,
        }
    }
}
