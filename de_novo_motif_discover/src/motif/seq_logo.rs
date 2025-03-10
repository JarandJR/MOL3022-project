use std::process::Command;

use super::pwm::PositionWeightMatrix;

pub fn generate_sequence_logo(pwm: PositionWeightMatrix) {
    let pwm_json = pwm.as_json();
    let output = Command::new("python")
        .arg("src/motif/plot_logo.py")
        .arg(pwm_json)
        .output()
        .expect("Failed to execute Python script");

    println!("Generating Sequence logo..");
    if output.status.success() {
        println!("Sequence logo generated successfully!");
    } else {
        eprintln!("Error: {}", String::from_utf8_lossy(&output.stderr));
    }
}
