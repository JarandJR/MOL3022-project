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
}

impl Default for PositionWeightMatrix {
    fn default() -> Self {
        // Testdata
        let pwm = vec![
            vec![0.2, 0.1, 0.4, 0.3, 0.2, 0.1, 0.4, 0.3], // A
            vec![0.1, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2], // C
            vec![0.1, 0.3, 0.2, 0.4, 0.1, 0.3, 0.2, 0.4], // G
            vec![0.6, 0.2, 0.3, 0.1, 0.4, 0.2, 0.3, 0.1], // T
        ];
        Self { pwm }
    }
}
