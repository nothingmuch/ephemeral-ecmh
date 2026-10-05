# sage/kat.sage for the 14-byte families over p = 2^107 - 1: certificates
# (src/curvegen/select107.rs) plus hash-to-curve and sum vectors on each certified
# curve, and each family's per-window rho cost on stderr.
#
#   nix develop .#math -c sage sage/kat107.sage > tests/common/kats107.rs
#   nix fmt
#
# The field and group arithmetic is independent of the Rust code; the candidate
# derivation, the byte framing (src/hash.rs, src/curvegen/select107.rs, the
# 108-bit encodings of src/field/fp107.rs) and the witness-selection rules match
# Rust by design. Rejection witnesses are the Rust prover's hashed points
# (sage/kat_common.sage), so a run reproduces every byte.
load("sage/kat_common.sage")

p = 2^107 - 1
fp = PrimeCodec(GF(p), 107)
TAG_FP = b"ephemeral-ecmh/curve/fp107"
TAG_WEIER = b"ephemeral-ecmh/curve/weier107"

# --- emit ------------------------------------------------------------------
out = sys.stdout
write_struct(
    out,
    "kat107.sage",
    "tests/common/kats107.rs",
    struct_doc=["/// As `kats::CertKat`, with 108-bit encodings (14 little-endian bytes)."],
)
out.write("\n")


def rho_bits(r):
    """log2 of the expected rho cost with the negation map, sqrt(pi r / 4)."""
    return float(log(sqrt(pi * r / 4), 2))


def emit(name, prove, cofactor):
    out.write("pub const %s: &[CertKat] = &[\n" % name)
    for seed in SEEDS:
        E, i, r, rej = prove(seed)
        assert E.order() == cofactor * r
        Ps = [fp.hash(E, seed, bytes([k])) for k in range(NMSG)]
        write_cert(out, seed, i, r, rej, [fp.encode(P) for P in Ps], fp.encode(sum(Ps)), 28)
        print(
            "%s: seed %s.. index %d, %d rejections, log2 r = %.2f, rho 2^%.2f"
            % (name, seed[:4].hex(), i, len(rej), float(log(r, 2)), rho_bits(r)),
            file=sys.stderr,
        )
    out.write("];\n\n")


emit("FP107_CERTS", lambda seed: edwards_prove(fp, TAG_FP, seed), 4)
emit("WEIER107_CERTS", lambda seed: weier_prove(fp, TAG_WEIER, seed), 1)
