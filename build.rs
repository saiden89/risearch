//! Declares the `kani` cfg so proofs compile without unexpected-cfg warnings.

fn main() {
    println!("cargo::rustc-check-cfg=cfg(kani)");
}
