use super::Sequence;

pub struct FastaSeq;

impl FastaSeq {
    pub fn sequences(contents: String) -> Vec<Sequence> {
        contents
            .lines()
            .skip(1)
            .step_by(2)
            .map(|s| s.to_owned())
            .collect()
    }
}
