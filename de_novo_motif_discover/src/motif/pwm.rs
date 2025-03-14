use std::fmt::Debug;

use crate::dna::nucleotide::Nucleotide;

pub struct PositionWeightMatrix {
    pwm: Vec<Vec<f64>>,
}

impl PositionWeightMatrix {
    pub fn new(pwm: Vec<Vec<f64>>) -> Self {
        Self { pwm }
    }

    pub fn as_json(&self) -> String {
        serde_json::to_string(&self.pwm).unwrap()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Vec<f64>> {
        self.pwm.iter()
    }
}

impl From<Vec<Vec<f64>>> for PositionWeightMatrix {
    fn from(value: Vec<Vec<f64>>) -> Self {
        PositionWeightMatrix::new(value)
    }
}

impl Debug for PositionWeightMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\n")?;
        self.pwm
            .iter()
            .enumerate()
            .map(|(nuc, pwm)| write!(f, "{}: {:?}\n", Nucleotide::from(nuc), pwm))
            .collect()
    }
}
