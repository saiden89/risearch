//! Code-generates the built-in DSM tables from `data/dsm/manifest.toml`.

use serde::Deserialize;
use std::{
    env, fs,
    path::{Path, PathBuf},
};

#[derive(Deserialize)]
struct Manifest {
    dsm: Vec<DsmEntry>,
}

#[derive(Deserialize)]
struct DsmEntry {
    id: String,
    orientation: String,
    invalid_transition_kcal: f64,
    temperatures: Vec<TemperatureEntry>,
}

#[derive(Deserialize)]
struct TemperatureEntry {
    temperature: i32,
    file: PathBuf,
    initiation_kcal: f64,
}

fn main() {
    println!("cargo::rustc-check-cfg=cfg(kani)");
    println!("cargo:rerun-if-changed=data/dsm/");
    let manifest_path = Path::new("data/dsm/manifest.toml");
    let manifest_content = fs::read_to_string(manifest_path).expect("Failed to read manifest");
    let manifest: Manifest = toml::from_str(&manifest_content).expect("Failed to parse manifest");

    let mut builtin_tables = String::from("static BUILTIN_TABLES: &[Dsm] = &[\n");
    let mut builtin_names = String::from("static BUILTIN_NAMES: &[&str] = &[\n");

    for entry in &manifest.dsm {
        builtin_names.push_str(&format!("    \"{}\",\n", entry.id));

        let orient = match entry.orientation.as_str() {
            "identity" => "Orientation::Identity",
            "reverse-swap" => "Orientation::ReverseSwap",
            _ => panic!("unknown orientation"),
        };

        for temp in &entry.temperatures {
            let tsv = fs::canonicalize(Path::new("data/dsm").join(&temp.file))
                .expect("Failed to canonicalize TSV path");
            let tsv_str = tsv.to_str().expect("TSV path is not UTF-8");

            builtin_tables.push_str("    Dsm {\n");
            builtin_tables.push_str(&format!("        id: \"{}\",\n", entry.id));
            builtin_tables.push_str(&format!("        temperature: {},\n", temp.temperature));
            builtin_tables.push_str(&format!(
                "        initiation: {:?},\n",
                temp.initiation_kcal
            ));
            builtin_tables.push_str(&format!(
                "        invalid_transition: {:?},\n",
                entry.invalid_transition_kcal
            ));
            builtin_tables.push_str(&format!("        orientation: {},\n", orient));
            builtin_tables.push_str(&format!(
                "        tsv_content: include_str!(r\"{}\"),\n",
                tsv_str
            ));
            builtin_tables.push_str("    },\n");
        }
    }

    builtin_tables.push_str("];\n");
    builtin_names.push_str("];\n");

    let out_dir = env::var("OUT_DIR").expect("OUT_DIR not set");
    let dest_path = Path::new(&out_dir).join("generated_canonical_tables.rs");
    fs::write(dest_path, format!("{}\n{}", builtin_tables, builtin_names))
        .expect("Failed to write generated file");
}
