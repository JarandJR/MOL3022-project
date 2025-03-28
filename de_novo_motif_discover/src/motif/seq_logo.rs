use std::process::Command;

use crate::methods::SupportedMethods;

use super::pwm::PositionWeightMatrix;

pub fn generate_sequence_logo(pwm: &PositionWeightMatrix, method: &SupportedMethods) {
    println!("Generating Sequence logo..");
    let pwm_json = pwm.as_json();
    let output = Command::new("bash")  // Changed from "sh" to "bash"
        .arg("-c")
        .arg(format!(". .env/bin/activate && python src/motif/plot_logo.py '{}' '{}'", 
                     pwm_json, 
                     (format!("{:?}", method)).replace(' ', "_").to_lowercase()))
        .output()
        .expect("Failed to execute Python script");

    if output.status.success() {
        println!("Sequence logo generated successfully!");
    } else {
        eprintln!("Error: {}", String::from_utf8_lossy(&output.stderr));
    }
}
