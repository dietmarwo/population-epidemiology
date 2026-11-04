fn main() {
    let rustc = std::env::var_os("RUSTC").expect("Cargo supplies RUSTC");
    let output = std::process::Command::new(rustc)
        .arg("--version")
        .output()
        .expect("rustc version");
    assert!(output.status.success());
    println!(
        "cargo:rustc-env=CRN_RUSTC_VERSION={}",
        String::from_utf8(output.stdout)
            .expect("UTF-8 rustc version")
            .trim()
    );
    println!("cargo:rerun-if-changed=../../Cargo.lock");
    println!("cargo:rerun-if-changed=build.rs");
}
