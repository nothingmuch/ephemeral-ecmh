# Independent quadratic-field selector/certificate and codec reference.
#
#   nix develop .#math -c python sage/kat_fp2.sage --self-test
#   nix develop .#math -c python sage/kat_fp2.sage > tests/common/kats_fp2.rs
#
# Python integer-pair arithmetic, generic Tonelli-Shanks, and complete affine
# laws are independent of Rust. PARI supplies ellcard, factor, and proved
# primality: counts share Rust's backend and are binding checks, not an
# independent point-counting algorithm. No Rust library is called here.
import argparse
import atexit
import hashlib
import sys
import subprocess


TAG_CERT = b"ephemeral-ecmh/cert-point"
TAG_WITNESS = b"ephemeral-ecmh/prove/witness"
SEEDS = [bytes.fromhex("000000000019d6689c085ae165831e934ff763ae46a2a6c172b3f1b60a8ce26f")]
INVALID = object()


def salted(tag, seed, msg, counter):
    t = hashlib.sha256(tag).digest()
    return hashlib.sha256(t + t + seed + bytes(32) + msg + counter.to_bytes(4, "little")).digest()


def halves(digest):
    return [int.from_bytes(digest[:16], "little"), int.from_bytes(digest[16:], "little")]


class Quadratic:
    """Canonical a+b*i, i^2=nu, with the specified coefficient packing."""
    def __init__(self, p, nu, width):
        self.p, self.nu, self.width = p, nu % p, width
        self.q, self.bits = p * p, 2 * width
        self.mask = (1 << self.bits) - 1
        self._nonsquare = None

    def __call__(self, a=0, b=0):
        if isinstance(a, Element):
            assert a.f is self and b == 0
            return a
        return Element(self, a % self.p, b % self.p)

    def unpack(self, value):
        if not 0 <= value <= self.mask:
            return None
        a, b = value & ((1 << self.width) - 1), value >> self.width
        return self(a, b) if a < self.p and b < self.p else None

    def reduce(self, value):
        value &= self.mask
        return self(value & ((1 << self.width) - 1), value >> self.width)

    def nonsquare(self):
        if self._nonsquare is None:
            # All base-field elements are squares in this quadratic field.
            for a in range(self.p):
                z = self(a, 1)
                if not z.is_square():
                    self._nonsquare = z
                    break
        return self._nonsquare


class Element:
    __slots__ = ("f", "a", "b")

    def __init__(self, f, a, b):
        self.f, self.a, self.b = f, a, b

    def __repr__(self):
        return "(%d+%d*i)" % (self.a, self.b)

    def __hash__(self):
        return hash((id(self.f), self.a, self.b))

    def __eq__(self, other):
        if isinstance(other, int):
            other = self.f(other)
        return isinstance(other, Element) and self.f is other.f and (self.a, self.b) == (other.a, other.b)

    def __neg__(self):
        return self.f(-self.a, -self.b)

    def __add__(self, other):
        other = self.f(other)
        return self.f(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __sub__(self, other):
        return self + -self.f(other)

    def __rsub__(self, other):
        return self.f(other) - self

    def __mul__(self, other):
        other = self.f(other)
        return self.f(self.a * other.a + self.f.nu * self.b * other.b,
                      self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def inverse(self):
        norm = (self.a * self.a - self.f.nu * self.b * self.b) % self.f.p
        assert norm != 0, "inverse of zero"
        inv = pow(norm, -1, self.f.p)
        return self.f(self.a * inv, -self.b * inv)

    def __truediv__(self, other):
        return self * self.f(other).inverse()

    def __pow__(self, exponent):
        if exponent < 0:
            return self.inverse() ** -exponent
        result, base = self.f(1), self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def pack(self):
        return self.a | (self.b << self.f.width)

    def sign(self):
        return (self.a if self.a else self.b) & 1

    def is_square(self):
        return self == 0 or self ** ((self.f.q - 1) // 2) == 1

    def sqrt(self):
        """Generic Tonelli-Shanks over Fp^2, unlike Rust's complex method."""
        if self == 0:
            return self
        if not self.is_square():
            return None
        odd, shift = self.f.q - 1, 0
        while odd % 2 == 0:
            odd //= 2
            shift += 1
        c = self.f.nonsquare() ** odd
        x, t = self ** ((odd + 1) // 2), self ** odd
        while t != 1:
            i, tt = 0, t
            while tt != 1:
                tt = tt * tt
                i += 1
            assert 0 < i < shift
            b = c ** (1 << (shift - i - 1))
            x, t, c, shift = x * b, t * b * b, b * b, i
        assert x * x == self
        return -x if x.sign() else x


class Curve:
    """Twisted raw Edwards points, or full cubic-model affine points."""
    def __init__(self, family, parameter):
        self.family, self.f, self.parameter = family, family.field, parameter
        self.quotient = family.model == "twisted"
        if family.model == "weier":
            self.a2, self.a4, self.a6 = self.f(0), self.f(-3), parameter
        else:
            a = -1 if self.quotient else 1
            self.a2 = (a + parameter) / 2
            self.a4 = ((a - parameter) / 4) ** 2
            self.a6 = self.f(0)
        self.zero = (self.f(0), self.f(1)) if self.quotient else None

    def j(self):
        if self.family.model == "weier":
            return self.f(6912) / (4 - self.parameter ** 2)
        a, b = self.a2, self.a4
        return 256 * (a * a - 3 * b) ** 3 / (b * b * (a * a - 4 * b))

    def identity(self, point):
        return point[0] == 0 if self.quotient else point is None

    def on_curve(self, point):
        if not self.quotient and point is None:
            return True
        x, y = point
        if self.quotient:
            return y * y - x * x == 1 + self.parameter * x * x * y * y
        return y * y == self.rhs(x)

    def rhs(self, x):
        return x * (x * (x + self.a2) + self.a4) + self.a6

    def negate(self, point):
        if point is None:
            return None
        x, y = point
        return (-x, y) if self.quotient else (x, -y)

    def add(self, left, right):
        if self.quotient:
            u, v = left
            x, y = right
            t = self.parameter * u * v * x * y
            return ((u * y + v * x) / (1 + t), (v * y + u * x) / (1 - t))
        if left is None:
            return right
        if right is None:
            return left
        x, y = left
        xx, yy = right
        if x == xx:
            if y == -yy:
                return None
            slope = (3 * x * x + 2 * self.a2 * x + self.a4) / (2 * y)
        else:
            slope = (yy - y) / (xx - x)
        xxx = slope * slope - self.a2 - x - xx
        return xxx, slope * (x - xxx) - y

    def mul(self, n, point):
        result = self.zero
        while n:
            if n & 1:
                result = self.add(result, point)
            point = self.add(point, point)
            n >>= 1
        return result

    def encode(self, point):
        if self.quotient:
            u, v = point
            flip = u.sign() if v == 0 else v.sign()
            return (-u if flip else u).pack()
        if point is None:
            return self.f.mask
        x, y = point
        return x.pack() | (y.sign() << self.f.bits)

    def decode(self, code):
        if self.quotient:
            u = self.f.unpack(code)
            if u is None:
                return INVALID
            den = 1 - self.parameter * u * u
            assert den != 0
            v = ((1 + u * u) / den).sqrt()
            if v is None or (v == 0 and u.sign()):
                return INVALID
            return u, v
        if not 0 <= code < (1 << (self.f.bits + 1)):
            return INVALID
        packed, sign = code & self.f.mask, code >> self.f.bits
        if packed == self.f.mask:
            return INVALID if sign else None
        x = self.f.unpack(packed)
        if x is None:
            return INVALID
        y = self.rhs(x).sqrt()
        if y is None or (y == 0 and sign):
            return INVALID
        return x, -y if y.sign() != sign else y

    def hash_point(self, seed, msg, tag=b"kat"):
        bits = self.f.bits if self.quotient else self.f.bits + 1
        mask = (1 << bits) - 1
        for counter in range(2**32):
            for half in halves(salted(tag, seed, msg, counter)):
                code = half & mask
                if self.f.unpack(code & self.f.mask) is None:
                    continue
                point = self.decode(code)
                if point is not INVALID:
                    assert self.on_curve(point) and self.encode(point) == code
                    return point
        raise RuntimeError("point hash counter exhausted")


class Family:
    def __init__(self, name, tag, p, nu, width, model):
        self.name, self.tag, self.model = name, tag, model
        self.field = Quadratic(p, nu, width)
        self.h = 1 if model == "weier" else 4
        self.lo, self.hi = (p - 1)**2, (p + 1)**2
        self.r_lo = self.lo // self.h
        assert self.r_lo > self.hi - self.lo

    def construct(self, parameter, exclude_subfield=True):
        if self.model == "weier":
            if parameter in (0, 2, -2):
                return None
        elif parameter.is_square():
            return None
        curve = Curve(self, parameter)
        if exclude_subfield and curve.j().b == 0:
            return None
        return curve

    def candidate(self, seed, index):
        raw = halves(salted(self.tag, seed, b"", index))[0]
        return self.construct(self.field.reduce(raw))


FAMILIES = [
    Family("TWISTED61X2", b"ephemeral-ecmh/curve/twisted61x2", 2**61-1, -1, 61, "twisted"),
    Family("TWISTED64X2", b"ephemeral-ecmh/curve/twisted64x2", 2**64-59, 2, 64, "twisted"),
    Family("TWISTED_GOLDILOCKS2", b"ephemeral-ecmh/curve/twisted-goldilocks2", 2**64-2**32+1, 7, 64, "twisted"),
    Family("EDWARDS61X2", b"ephemeral-ecmh/curve/edwards61x2", 2**61-1, -1, 61, "edwards"),
    Family("WEIER61X2", b"ephemeral-ecmh/curve/weier61x2", 2**61-1, -1, 61, "weier"),
]


class Backend:
    """A single persistent PARI/GP process; no Sage/cypari2 dependency."""
    def __init__(self, executable="gp"):
        self.process = subprocess.Popen([executable, "-q", "-f", "-s", "64000000",
                                         "-D", "parisizemax=4000000000"], stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE, text=True, bufsize=1)
        atexit.register(self.close)
        version = self.evaluate("version()")
        print("BACKEND PARI %s executable=%s" % (version, executable), file=sys.stderr, flush=True)
        self.generators = {}
        sieve = bytearray([1]) * 65536
        sieve[:2] = b"\x00\x00"
        for p in range(2, 256):
            if sieve[p]:
                sieve[p*p::p] = bytes(len(sieve[p*p::p]))
        self.small_primes = [p for p in range(3, 65536, 2) if sieve[p]]

    def evaluate(self, expression):
        # The second marker also follows a PARI exception, avoiding a blocked
        # reader if one of these generated, integer-only expressions fails.
        request = 'iferr(fp2_result=(%s);print("FP2_RESULT:",fp2_result),fp2_error,print("FP2_ERROR:",Str(fp2_error)));print("FP2_DONE");\n' % expression
        self.process.stdin.write(request)
        self.process.stdin.flush()
        result, error = None, None
        for line in self.process.stdout:
            line = line.strip()
            if line.startswith("FP2_RESULT:"):
                result = line[len("FP2_RESULT:"):]
            elif line.startswith("FP2_ERROR:"):
                error = line[len("FP2_ERROR:"):]
            elif line == "FP2_DONE":
                if error is not None:
                    raise RuntimeError("PARI: " + error)
                assert result is not None, "PARI returned no result"
                return result
        raise RuntimeError("PARI process exited unexpectedly")

    def order(self, curve):
        f = curve.f
        if f not in self.generators:
            name = "fp2_generator_%d" % len(self.generators)
            expression = "%s=ffgen(Mod(1,%d)*(x^2-%d))" % (name, f.p, f.nu)
            assert self.evaluate("type(%s)" % expression) == "t_FFELT"
            self.generators[f] = name
        name = self.generators[f]
        value = lambda x: "(%d+%d*%s)" % (x.a, x.b, name)
        coefficients = "[0,%s,0,%s,%s]" % (value(curve.a2), value(curve.a4), value(curve.a6))
        return int(self.evaluate("ellcard(ellinit(%s,%s))" % (coefficients, name)))

    def isprime(self, n):
        return bool(int(self.evaluate("isprime(%d)" % n)))

    def znorder(self, q, r):
        """The multiplicative order of q mod r: the embedding degree of r over F_q."""
        return int(self.evaluate("znorder(Mod(%d, %d))" % (q, r)))

    def factor(self, n):
        for p in self.small_primes:
            if n % p == 0:
                return p
        return int(self.evaluate("factor(%d)[1,1]" % n))

    def close(self):
        if self.process.poll() is None:
            try:
                self.process.stdin.write("quit;\n")
                self.process.stdin.flush()
                self.process.wait(timeout=5)
            except (BrokenPipeError, subprocess.TimeoutExpired):
                self.process.kill()
                self.process.wait()


def witness(curve, order, label, seed, index):
    """Prime-power extraction in E, or in the quotient of order #E/2."""
    base = 2 if label in (2, 8) else label
    target = 4 if curve.quotient and label == 8 else label
    size = order // 2 if curve.quotient else order
    reduced = size
    while reduced % base == 0:
        reduced //= base
    for attempt in range(2**32):
        msg = index.to_bytes(4, "little") + attempt.to_bytes(4, "little")
        point = curve.mul(reduced, curve.hash_point(seed, msg, TAG_WITNESS))
        if curve.identity(point):
            continue
        while not curve.identity(curve.mul(target, point)):
            point = curve.mul(base, point)
        if curve.identity(curve.mul(target // base, point)):
            continue
        code = curve.encode(point)
        decoded = curve.decode(code)
        assert decoded is not INVALID and curve.on_curve(decoded)
        assert curve.encode(decoded) == code
        assert not curve.identity(curve.mul(target // base, decoded))
        assert curve.identity(curve.mul(target, decoded))
        return label, code
    raise RuntimeError("witness counter exhausted")


def acceptance(curve, seed, index, r):
    point = curve.hash_point(seed, index.to_bytes(4, "little"), TAG_CERT)
    point = curve.mul(2 if curve.quotient else curve.family.h, point)
    assert not curve.identity(point), "certificate hash landed in the cofactor subgroup"
    assert curve.identity(curve.mul(r, point))


def prove(family, seed, backend):
    rejections, orders, factors = [], [], []
    for index in range(2**32):
        curve = family.candidate(seed, index)
        if curve is None:
            continue
        order = backend.order(curve)
        assert family.lo <= order <= family.hi
        assert order % family.h == 0
        orders.append((curve.parameter.pack(), order))
        print("COUNT %s parameter=%d order=%d" % (family.name, curve.parameter.pack(), order), file=sys.stderr, flush=True)
        r = order // family.h
        special = (2 if order % 2 == 0 else None) if family.h == 1 else (8 if order % 8 == 0 else None)
        if special is None and backend.isprime(r):
            assert family.r_lo <= r <= family.hi // family.h
            assert backend.znorder(family.field.q, r) > 2**20, "embedding degree at most 2^20"
            acceptance(curve, seed, index, r)
            return curve, index, r, rejections, orders, factors
        if special is None:
            label = backend.factor(r)
            assert 3 <= label < family.r_lo and r % label == 0 and backend.isprime(label)
            factors.append((r, label))
        else:
            label = special
        rejections.append(witness(curve, order, label, seed, index))
        print("%s seed=%s index=%d rejected by %d" % (family.name, seed.hex()[:12], index, label), file=sys.stderr, flush=True)
    raise RuntimeError("candidate counter exhausted")


def verify_rejections(family, seed, index, rejections, candidate=None):
    """Structural skips consume no entry; an injected schedule tests this."""
    candidate = family.candidate if candidate is None else candidate
    position = 0
    for k in range(index):
        curve = candidate(seed, k)
        if curve is None:
            continue
        assert position < len(rejections)
        label, code = rejections[position]
        position += 1
        point = curve.decode(code)
        assert point is not INVALID and not curve.identity(point)
        if label == 8 and family.h == 4:
            target = 4 if curve.quotient else 8
            assert not curve.identity(curve.mul(target // 2, point))
            assert curve.identity(curve.mul(target, point))
        elif label == 2 and family.h == 1:
            assert point[1] == 0
        else:
            assert label >= 3 and label % 2 and label < family.r_lo
            assert curve.identity(curve.mul(label, point))
    assert position == len(rejections)


def verify_reference(family, seed, index, r, rejections, backend):
    verify_rejections(family, seed, index, rejections)
    assert family.lo <= family.h * r <= family.hi and r > family.hi - family.lo
    assert backend.znorder(family.field.q, r) > 2**20
    curve = family.candidate(seed, index)
    assert curve is not None
    acceptance(curve, seed, index, r)


def self_test():
    """Boundary cases and tiny-field enumeration need no PARI backend."""
    for family in FAMILIES:
        f = family.field
        assert f(0, 1)**2 == f(f.nu)
        assert f.unpack(f.mask) is None
        assert f.unpack(f.p) is None and f.unpack(f.p << f.width) is None
        assert f.unpack((1 << f.bits)) is None
        assert f.reduce(f.p) == 0 and f.reduce(f.p << f.width) == 0
        for x in (f(0), f(1), f(0, 1), f(f.p-1, f.p-2)):
            assert f.unpack(x.pack()) == x
            assert x == 0 or x.sign() != (-x).sign()
            assert (x*x).sqrt()**2 == x*x
        if family.model == "weier":
            outside = f(0, 1)
            curve = family.construct(outside, exclude_subfield=False)
            assert outside.b != 0 and curve.j().b == 0 and curve.j()**f.p == curve.j()
            assert family.construct(outside) is None
            assert family.construct(f(1)) is None
        curve = next(c for k in range(100) if (c := family.candidate(bytes(32), k)) is not None)
        assert curve.j().b != 0 and curve.j()**f.p != curve.j()
        assert curve.decode(1 << (f.bits + (0 if curve.quotient else 1))) is INVALID
        assert curve.decode(f.p) is INVALID
        assert curve.decode(f.p << f.width) is INVALID
        if curve.quotient:
            assert curve.encode((f(0), f(1))) == curve.encode((f(0), f(-1))) == 0
            i = f(-1).sqrt()
            assert i is not None
            for point in ((i, f(0)), (-i, f(0))):
                code = curve.encode(point)
                assert curve.decode(code) is not INVALID
                assert curve.identity(curve.mul(2, point)) and not curve.identity(point)
            assert curve.decode((-i if i.sign() == 0 else i).pack()) is INVALID
        else:
            assert curve.decode(f.mask) is None
            assert curve.decode(f.mask | (1 << f.bits)) is INVALID
            if family.model == "edwards":
                assert curve.decode(0) == (f(0), f(0))
                assert curve.decode(1 << f.bits) is INVALID
        points = [curve.hash_point(bytes(32), bytes([k])) for k in range(3)]
        for point in points:
            assert curve.on_curve(point)
            assert curve.identity(curve.add(point, curve.negate(point)))
            assert curve.encode(curve.decode(curve.encode(point))) == curve.encode(point)
    # A real non-base coefficient with subfield j is skipped before a real
    # order-two rejection. Missing and extra entries must both fail.
    family = FAMILIES[-1]
    f = family.field
    skipped = family.construct(f(0, 1))
    assert skipped is None
    x = f(1, 1)
    curve = family.construct(3*x - x**3)
    assert curve is not None and curve.rhs(x) == 0
    rejection = (2, curve.encode((x, f(0))))
    assert curve.decode(rejection[1] | (1 << f.bits)) is INVALID
    schedule = lambda seed, index: [skipped, curve][index]
    verify_rejections(family, bytes(32), 2, [rejection], schedule)
    for entries in ([], [rejection, rejection]):
        try:
            verify_rejections(family, bytes(32), 2, entries, schedule)
        except AssertionError:
            pass
        else:
            raise AssertionError("skip traversal accepted a missing or extra entry")
    # Exhaust the quadratic field F_3[i]/(i^2+1) and every complete twisted
    # parameter. Check quotient codec cosets, closure, and associativity.
    tiny = Family.__new__(Family)
    tiny.name, tiny.model, tiny.field = "TINY", "twisted", Quadratic(3, -1, 2)
    f = tiny.field
    elements = [f(a,b) for a in range(3) for b in range(3)]
    for d in elements:
        if d.is_square():
            continue
        curve = Curve(tiny, d)
        points = [(u,v) for u in elements for v in elements if curve.on_curve((u,v))]
        encodings = {curve.encode(point) for point in points}
        assert len(encodings)*2 == len(points)
        for point in points:
            decoded = curve.decode(curve.encode(point))
            assert decoded in (point, (-point[0], -point[1]))
            for other in points:
                assert curve.on_curve(curve.add(point, other))
                decoded_other = curve.decode(curve.encode(other))
                assert curve.encode(curve.add(point, other)) == curve.encode(curve.add(decoded, decoded_other))
                assert curve.encode(curve.add(point, curve.negate(other))) == curve.encode(curve.add(decoded, curve.negate(decoded_other)))
                for third in points:
                    assert curve.add(curve.add(point, other), third) == curve.add(point, curve.add(other, third))
        for code in range(1 << f.bits):
            assert (curve.decode(code) is not INVALID) == (code in encodings)
    print("PASS: exact bases, packing/sign boundaries, structural j exclusions, quotient kernels, affine laws", file=sys.stderr)


def emit_header():
    print("// @generated by sage/kat_fp2.sage; do not edit. Regenerate from the repository root:")
    print("//   nix develop .#math -c python sage/kat_fp2.sage > tests/common/kats_fp2.rs && nix fmt")
    print("// Independent Python field/affine reference; point counts and proved primes use PARI.")
    print("// PARI is also Rust's count backend, so counts are binding checks, not independent counts.")
    print("#![allow(dead_code)]\n")
    print("pub struct CertKat {\n    pub seed: [u8; 32],\n    pub index: u32,\n    pub r: u128,")
    print("    pub rejections: &'static [(u128, u128)],")
    print('    /// Salted(b"kat", seed), messages [i] for i in 0..8.\n    pub hashes: [u128; 8],')
    print("    pub sum: u128,\n    /// hash[0] - hash[1] + ... - hash[7].\n    pub signed_sum: u128,")
    print("    /// Canonical parameter and full #E, including the accepted candidate.\n    pub orders: &'static [(u128, u128)],")
    print("    /// Composite odd m and its smallest prime factor.\n    pub factors: &'static [(u128, u128)],\n}\n")


def emit(families, executable):
    backend = Backend(executable)
    emit_header()
    for family in families:
        print("pub const %s_CERTS: &[CertKat] = &[" % family.name)
        coverage = set()
        seeds = list(SEEDS)
        for number, seed in enumerate(seeds):
            curve, index, r, rejected, orders, factors = prove(family, seed, backend)
            verify_reference(family, seed, index, r, rejected, backend)
            coverage.update(label for label, _ in rejected)
            points = [curve.hash_point(seed, bytes([k])) for k in range(8)]
            total, signed = curve.zero, curve.zero
            for k, point in enumerate(points):
                total = curve.add(total, point)
                signed = curve.add(signed, curve.negate(point) if k % 2 else point)
            print("    CertKat {\n        seed: %s,\n        index: %d,\n        r: 0x%x," % (list(seed), index, r))
            print("        rejections: &[")
            for label, code in rejected:
                print("            (0x%x, 0x%032x)," % (label, code))
            print("        ],\n        hashes: [")
            for point in points:
                print("            0x%032x," % curve.encode(point))
            print("        ],\n        sum: 0x%032x,\n        signed_sum: 0x%032x," % (curve.encode(total), curve.encode(signed)))
            for name, pairs in (("orders", orders), ("factors", factors)):
                print("        %s: &[" % name)
                for a, b in pairs:
                    print("            (0x%x, 0x%x)," % (a, b))
                print("        ],")
            print("    },", flush=True)
            print("PASS %s: index=%d, rejections=%d, r proven prime" % (family.name, index, len(rejected)), file=sys.stderr, flush=True)
            complete = (2 if family.h == 1 else 8) in coverage and any(label >= 3 and label % 2 for label in coverage)
            if not complete and number == 0:
                seeds.append(hashlib.sha256(b"ephemeral-ecmh/kat/fp2/" + family.name.encode() + b"/1").digest())
        assert (2 if family.h == 1 else 8) in coverage, "need another seed for even rejection coverage"
        assert any(label >= 3 and label % 2 for label in coverage), "need another seed for odd rejection coverage"
        print("];\n", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--gp", default="gp", help="PARI/GP executable for counts and proved primality")
    parser.add_argument("--family", choices=[f.name for f in FAMILIES], action="append")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        emit([f for f in FAMILIES if args.family is None or f.name in args.family], args.gp)
