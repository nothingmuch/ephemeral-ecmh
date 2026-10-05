"""Bounded search for complete Edwards curves over pseudo-Mersenne primes.

Run with ``sage -python sage/edwards_search.py``. PARI proves primality,
counts points and factors orders; the script records every counted candidate.
The default domain is p = 2**bits - c for bits in {127, 128}, 3 <= c <= 1000,
and the first proven prime with p = 5 (mod 8) at each size. Eight seeded
nonsquares are sampled from the full field. Distinct nonsquares with d or 2d in
{1, -1, ..., 16, -16} form a separate small-constant comparison.

For a = -1, p = 1 (mod 4) and nonsquare d give complete Edwards addition.
The birational Montgomery model B*v**2 = u**3 + A*u**2 + u has
A = 2*(a+d)/(a-d), B = 4/(a-d). Scaling x = u/B and y = v/B gives the
Weierstrass coefficients used below. See HWCD, ePrint 2008/522, section 2.

The Decaf width is a capacity observation, not an implemented codec or a
performance result. Candidate screening is not a complete security analysis.
"""

import argparse
import hashlib
import json
import os
from functools import lru_cache

from cypari2 import Pari


@lru_cache(maxsize=1)
def pari():
    engine = Pari()
    engine.allocatemem(64 << 20, silent=True)
    if os.environ.get("GP_DATA_DIR"):
        engine.default("datadir", os.environ["GP_DATA_DIR"])
    return engine


def prime_fields(bits, max_c, limit):
    """Yield proven primes of exactly ``bits`` bits in increasing c order."""
    if bits < 3 or max_c < 1 or limit < 1:
        raise ValueError("bits >= 3, max_c >= 1 and limit >= 1 are required")
    found = 0
    stop = min(max_c, (1 << (bits - 1)) - 1)
    for c in range(3, stop + 1, 8):
        p = (1 << bits) - c
        if pari()(p).isprime():
            yield p, c
            found += 1
            if found == limit:
                return


def small_d(p, bound):
    """Yield nonsquares in signed-magnitude order, with unique residues."""
    if p <= 3 or p % 8 != 5 or bound < 1:
        raise ValueError("p = 5 mod 8 and a positive d bound are required")
    seen = set()
    for magnitude in range(1, bound + 1):
        for d in (magnitude, -magnitude):
            residue = d % p
            if residue not in seen:
                seen.add(residue)
                if residue not in (0, p - 1) and pow(residue, (p - 1) // 2, p) == p - 1:
                    yield d


def cheap_parameters(p, bound):
    """Enumerate cheap d first, then cheap 2d, without repeated curves."""
    seen = set()
    for d in small_d(p, bound):
        seen.add(d % p)
        yield d, "d", d
    half = pow(2, -1, p)
    for magnitude in range(1, bound + 1):
        for delta in (magnitude, -magnitude):
            residue = delta * half % p
            if residue in seen:
                continue
            seen.add(residue)
            if residue not in (0, p - 1) and pow(residue, (p - 1) // 2, p) == p - 1:
                d = residue if residue <= p // 2 else residue - p
                yield d, "2d", delta


def random_d(p, seed, count):
    """Sample distinct nonsquares across the full field with a recorded seed.

    Rejection of out-of-range integers avoids modular-reduction bias. The
    field and counter are framed in the SHAKE input. The finite attempt bound
    makes exhaustion explicit rather than silently shortening the sample.
    """
    if p <= 3 or p % 8 != 5 or not pari()(p).isprime():
        raise ValueError("p must be a proven prime congruent to 5 mod 8")
    if count < 0 or count > (p - 1) // 2:
        raise ValueError("count exceeds the number of distinct nonsquares")
    width = (p.bit_length() + 7) // 8
    prefix = (
        b"ephemeral-ecmh/edwards-search/d/v1"
        + width.to_bytes(4, "little")
        + p.to_bytes(width, "little")
        + len(seed).to_bytes(4, "little")
        + seed
    )
    seen = set()
    for counter in range(128 + 16 * count):
        if len(seen) == count:
            return
        digest = hashlib.shake_256(prefix + counter.to_bytes(8, "little")).digest(width)
        d = int.from_bytes(digest, "little") & ((1 << p.bit_length()) - 1)
        if d >= p or d in seen or pow(d, (p - 1) // 2, p) != p - 1:
            continue
        seen.add(d)
        yield d
    if len(seen) != count:
        raise RuntimeError(
            "dense-parameter sampling exhausted its finite attempt bound"
        )


def candidate(p, d, *, min_prime_bits=112, embedding_bound=100):
    """Count and screen one complete a=-1 Edwards curve, retaining failures."""
    if p <= 3 or p % 8 != 5 or not pari()(p).isprime():
        raise ValueError("p must be a proven prime congruent to 5 mod 8")
    if d % p in (0, p - 1) or pow(d, (p - 1) // 2, p) != p - 1:
        raise ValueError("d must define a nonsingular complete Edwards curve")
    if min_prime_bits < 1 or embedding_bound < 0:
        raise ValueError("invalid screening bounds")
    a2 = (d - 1) * pow(2, -1, p) % p
    a4 = (d + 1) ** 2 * pow(16, -1, p) % p
    curve = pari().ellinit([0, a2, 0, a4, 0], p)
    order = int(curve.ellcard())
    matrix = pari()(order).factor()
    factors = [(int(matrix[i, 0]), int(matrix[i, 1])) for i in range(matrix.nrows())]
    if not all(pari()(prime).isprime() for prime, _ in factors):
        raise ArithmeticError("order factorization contains an unproven prime")
    if order % 4:
        raise ArithmeticError(
            "complete Edwards curve must have order divisible by four"
        )
    largest = max((prime for prime, _ in factors if prime != 2), default=1)
    degree = None
    power = 1
    if largest > 1:
        for k in range(1, embedding_bound + 1):
            power = power * p % largest
            if power == 1:
                degree = k
                break
    rejections = []
    if largest == 1 or largest.bit_length() < min_prime_bits:
        rejections.append("prime_factor_below_target")
    if order == p:
        rejections.append("anomalous")
    if degree is not None:
        rejections.append("embedding_degree_within_bound")
    bits = p.bit_length()
    return {
        "kind": "candidate",
        "p": str(p),
        "field_bits": bits,
        "c": str((1 << bits) - p),
        "a": -1,
        "d": str(d),
        "d_residue": str(d % p),
        "complete": True,
        "order": str(order),
        "factors": [[str(prime), exponent] for prime, exponent in factors],
        "primality_proven": True,
        "largest_odd_prime": str(largest),
        "large_prime_bits": largest.bit_length() if largest > 1 else 0,
        "cofactor_of_large_prime": str(order // largest),
        "two_adicity": (order & -order).bit_length() - 1,
        "embedding_degree": degree,
        "embedding_degree_bound": embedding_bound,
        "quotient_order": str(order // 4),
        "quotient_information_bits": (order // 4 - 1).bit_length(),
        "coordinate_sign_bits": bits + 1,
        "quotient_field_bits": bits,
        "quotient_field_bytes": (bits + 7) // 8,
        "quotient_codec": "Decaf capacity; codec not implemented here",
        "selected": not rejections,
        "rejections": rejections,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits", nargs="+", type=int, default=[127, 128])
    parser.add_argument("--max-c", type=int, default=1000)
    parser.add_argument("--primes-per-size", type=int, default=1)
    parser.add_argument("--d-bound", type=int, default=16)
    parser.add_argument("--random-curves", type=int, default=8)
    parser.add_argument("--seed-hex", default="00" * 32)
    parser.add_argument("--min-prime-bits", type=int, default=112)
    parser.add_argument("--embedding-bound", type=int, default=100)
    args = parser.parse_args()
    if (
        any(bits < 3 for bits in args.bits)
        or args.max_c < 1
        or args.primes_per_size < 1
        or args.d_bound < 0
        or args.random_curves < 0
        or args.min_prime_bits < 1
        or args.embedding_bound < 0
    ):
        parser.error("invalid search bounds")
    try:
        seed = bytes.fromhex(args.seed_hex)
    except ValueError:
        parser.error("seed-hex must contain hexadecimal bytes")
    print(
        json.dumps(
            {"kind": "parameters", "pari_version": str(pari().version()), **vars(args)}
        ),
        flush=True,
    )
    fields = counted = selected = 0
    for bits in dict.fromkeys(args.bits):
        for p, c in prime_fields(bits, args.max_c, args.primes_per_size):
            fields += 1
            print(
                json.dumps(
                    {
                        "kind": "field",
                        "p": str(p),
                        "bits": bits,
                        "c": c,
                        "primality_proven": True,
                    }
                ),
                flush=True,
            )
            controls = cheap_parameters(p, args.d_bound) if args.d_bound else []
            for d, coefficient, constant in controls:
                row = candidate(
                    p,
                    d,
                    min_prime_bits=args.min_prime_bits,
                    embedding_bound=args.embedding_bound,
                )
                row["cheap_coefficient"] = coefficient
                row["cheap_constant"] = constant
                row["experiment"] = "small_constant_control"
                row["namespace_note"] = (
                    "Small fixed palettes do not provide broad curve diversity"
                )
                print(json.dumps(row), flush=True)
                counted += 1
                selected += int(row["selected"])
            for index, d in enumerate(random_d(p, seed, args.random_curves)):
                row = candidate(
                    p,
                    d,
                    min_prime_bits=args.min_prime_bits,
                    embedding_bound=args.embedding_bound,
                )
                row["experiment"] = "seeded_dense_parameter"
                row["sample_index"] = index
                print(json.dumps(row), flush=True)
                counted += 1
                selected += int(row["selected"])
    print(
        json.dumps(
            {
                "kind": "summary",
                "fields": fields,
                "counted": counted,
                "selected": selected,
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
