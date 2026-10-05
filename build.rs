//! Links libpari for the `pari` feature. PARI_PREFIX (set by the dev shells
//! and the nix check) points at the install; without it the system linker
//! paths have to find it. Nothing is checked here, so type-checking builds
//! (clippy, doc) need no PARI at all.

fn main() {
    println!("cargo::rerun-if-env-changed=PARI_PREFIX");
    println!("cargo::rerun-if-env-changed=GP_DATA_DIR");
    if std::env::var_os("CARGO_FEATURE_PARI").is_none() {
        return;
    }
    if let Some(prefix) = std::env::var_os("PARI_PREFIX") {
        let lib = std::path::Path::new(&prefix).join("lib");
        println!("cargo::rustc-link-search=native={}", lib.display());
        println!("cargo::rustc-link-arg=-Wl,-rpath,{}", lib.display());
    }
    println!("cargo::rustc-link-lib=dylib=pari");
}
