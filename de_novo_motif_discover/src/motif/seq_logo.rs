use std::process::Command;

use crate::methods::SupportedMethods;

use super::pwm::PositionWeightMatrix;

pub fn generate_sequence_logo(pwm: &PositionWeightMatrix, method: &SupportedMethods) {
    let pwm_json = pwm.as_json();
    let output = Command::new("python")
        .arg("src/motif/plot_logo.py")
        .arg(pwm_json)
        .arg((format!("{:?}", method)).replace(' ', "_").to_lowercase())
        .output()
        .expect("Failed to execute Python script");

    println!("Generating Sequence logo..");
    if output.status.success() {
        println!("Sequence logo generated successfully!");
    } else {
        eprintln!("Error: {}", String::from_utf8_lossy(&output.stderr));
    }
}
