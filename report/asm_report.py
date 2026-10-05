"""The assembly probes (src/asm_probes.rs) out of rustc's --emit asm output,
one file per (target, CPU) configuration, summarised for comparison:

    asm_report.py OUT CONFIG=FILE.s [CONFIG=FILE.s ...]

CONFIG names the configuration, e.g. aarch64-linux-android/cortex-x4; its
target triple's architecture (the part before the first "-") picks the
mnemonic classes. Writes, under OUT:

    probes.tsv           one row per configuration and probe
    summary.md           closure instruction counts, probes by configuration,
                         and the instruction selection in each multiply
    CONFIG/PROBE.s       each probe's assembly, then its crate-local callees'

A probe's closure is the probe and every function it reaches by calls or
tail calls among the functions the file defines; calls out of the file
(memcpy, the allocator, panics) are listed instead. Counts are static:
each instruction counts once, regardless of loop iteration or branch
frequency. They measure code size, not execution time.

Symbols are demangled by the program named by --demangler (default
llvm-cxxfilt), which reads one name per line.
"""

import argparse
import re
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

PROBES = "ephemeral_ecmh::asm_probes::"

# a label opening a non-local symbol: ELF and Mach-O local labels start
# with .L / L / ltmp, which are branch targets, not functions
LABEL = re.compile(r"^([A-Za-z_$][\w.$]*):")
LOCAL = re.compile(r"^(\.L|L|ltmp)")
# identical functions the compiler merged: the later one an alias
ALIAS = re.compile(r"^([A-Za-z_$][\w.$]*) = ([A-Za-z_$][\w.$]*)$")

# mnemonic classes per architecture, on the mnemonic with any Apple
# arrangement suffix (pmull.1q, eor.16b) dropped
CLASSES = {
    "aarch64": {
        "mul": {
            "mul",
            "umulh",
            "smulh",
            "madd",
            "msub",
            "mneg",
            "umull",
            "umaddl",
            "umsubl",
            "smull",
        },
        "carry": {"adds", "adcs", "adc", "subs", "sbcs", "sbc", "ngcs", "negs"},
        "clmul": {"pmull", "pmull2"},
        "call": {"bl", "blr"},
    },
    "x86_64": {
        "mul": {"mul", "mulx", "imul"},
        "carry": {"adc", "sbb", "adcx", "adox", "setb", "setae"},
        "clmul": {"pclmulqdq", "vpclmulqdq"},
        "call": {"call"},
    },
}

# direct calls and jumps, and GOT entries: x86-64 ELF's sym@GOTPCREL,
# aarch64 ELF's :got:sym and Mach-O's sym@GOTPAGE
DIRECT = {"aarch64": {"bl", "b"}, "x86_64": {"call", "jmp"}}
GOT = {
    "aarch64": re.compile(r":got:([\w.$]+)|([\w.$]+)@GOTPAGE\b"),
    "x86_64": re.compile(r"([\w.$]+)@GOTPCREL"),
}

# indirect calls and jumps through a register name no symbol
X86_REGISTER = re.compile(r"^(r[abcd]x|r[sd]i|r[sb]p|r\d+)$")

COMMENT = {"aarch64": re.compile(r"\s*(//|;).*$"), "x86_64": re.compile(r"\s*#.*$")}


@dataclass
class Function:
    name: str
    lines: list[str] = field(default_factory=list)
    insns: list[tuple[str, str]] = field(default_factory=list)


def arch_of(config: str) -> str:
    arch = config.split("-", 1)[0]
    if arch not in CLASSES:
        raise ValueError(f"{config}: no mnemonic classes for {arch}")
    return arch


def parse(text: str, arch: str) -> tuple[dict[str, Function], dict[str, str]]:
    """The functions in an assembly listing by (raw) symbol, each from its
    label to its .cfi_endproc with its instructions as (mnemonic,
    operands), and the aliases among them, each to the function it names."""
    comment = COMMENT[arch]
    funcs: dict[str, Function] = {}
    aliases: dict[str, str] = {}
    label = None
    current = None
    for line in text.splitlines():
        m = ALIAS.match(line)
        if m:
            aliases[m.group(1)] = m.group(2)
            continue
        m = LABEL.match(line)
        if m:
            if not LOCAL.match(m.group(1)):
                label = m.group(1)
            if current is not None:
                current.lines.append(line)
            continue
        s = line.strip()
        if s == ".cfi_startproc" and label is not None:
            current = funcs[label] = Function(label, [f"{label}:"])
            continue
        if current is None:
            continue
        if s == ".cfi_endproc":
            current = None
            continue
        if not s or s.startswith("."):
            continue
        s = comment.sub("", s)
        if not s:
            continue
        current.lines.append(line)
        mnemonic, *operands = s.split(None, 1)
        current.insns.append((mnemonic, operands[0].strip() if operands else ""))
    for a in list(aliases):
        t, seen = aliases[a], {a}
        while t in aliases and t not in seen:
            seen.add(t)
            t = aliases[t]
        aliases[a] = t
    return funcs, aliases


def base(mnemonic: str) -> str:
    # b.eq is a conditional branch, not b with an arrangement
    return mnemonic if mnemonic.startswith("b.") else mnemonic.split(".")[0]


def targets(insns: list[tuple[str, str]], arch: str) -> list[str]:
    """The symbols an instruction list calls or tail-calls, directly or
    through an address loaded from the GOT (a call through a register, as
    for a function called more than once in position-independent code)."""
    out = []
    for mnemonic, operands in insns:
        b = base(mnemonic)
        syms = [g for m in GOT[arch].finditer(operands) for g in m.groups() if g]
        if b in DIRECT[arch] and not syms:
            # call sym@PLT, bl sym, b sym (a tail call)
            sym = operands.split("@")[0]
            if not (arch == "x86_64" and X86_REGISTER.match(sym)):
                syms = [sym]
        for sym in syms:
            if LABEL.match(sym + ":") and not LOCAL.match(sym) and sym not in out:
                out.append(sym)
    return out


def stats(insns: list[tuple[str, str]], arch: str) -> dict[str, int]:
    classes = CLASSES[arch]
    c = Counter()
    for mnemonic, operands in insns:
        b = base(mnemonic)
        for k, names in classes.items():
            if b in names:
                c[k] += 1
        if arch == "aarch64":
            if b.startswith("b.") or b in {"cbz", "cbnz", "tbz", "tbnz"}:
                c["branch"] += 1
            if b.startswith(("ld", "st")) and re.search(r"\[(sp|x29)\b", operands):
                c["stack"] += 1
        else:
            if b.startswith("j") and b != "jmp":
                c["branch"] += 1
            if b in {"push", "pop"} or re.search(r"\[(rsp|rbp)\b", operands):
                c["stack"] += 1
    c["insns"] = len(insns)
    return dict(c)


def demangle(names: list[str], demangler: str) -> dict[str, str]:
    # Mach-O prefixes every symbol with an underscore
    plain = [n[1:] if n.startswith(("__R", "__ZN")) else n for n in names]
    out = subprocess.run(
        [demangler],
        input="\n".join(plain) + "\n",
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    return dict(zip(names, out, strict=True))


def closure(
    funcs: dict[str, Function], aliases: dict[str, str], root: str, arch: str
) -> tuple[list[str], list[str]]:
    """The functions `root` reaches within the file, root first, and the
    symbols it calls outside it."""
    seen, external, todo = [root], [], [root]
    while todo:
        for t in targets(funcs[todo.pop()].insns, arch):
            t = aliases.get(t, t)
            if t in funcs:
                if t not in seen:
                    seen.append(t)
                    todo.append(t)
            elif t not in external:
                external.append(t)
    return seen, external


COLUMNS = [
    "insns",
    "closure",
    "callees",
    "branch",
    "stack",
    "mul",
    "carry",
    "clmul",
    "alias_of",
]


def probe_rows(config: str, text: str, demangler) -> tuple[list[dict], dict[str, str]]:
    """One row per probe in the listing, and each probe's assembly."""
    arch = arch_of(config)
    funcs, aliases = parse(text, arch)
    called = {t for f in funcs.values() for t in targets(f.insns, arch)}
    names = demangler(list(funcs) + sorted((called | set(aliases)) - set(funcs)))
    rows, listings = [], {}
    for raw in [*funcs, *aliases]:
        name = names[raw]
        if not name.startswith(PROBES):
            continue
        probe = name[len(PROBES) :]
        target = aliases.get(raw, raw)
        if target not in funcs:
            continue
        f = funcs[target]
        reach, external = closure(funcs, aliases, target, arch)
        own = stats(f.insns, arch)
        total = Counter()
        for r in reach:
            total.update(stats(funcs[r].insns, arch))
        rows.append(
            {
                "config": config,
                "probe": probe,
                "insns": own.get("insns", 0),
                "closure": total["insns"],
                "callees": len(reach) - 1,
                **{k: total[k] for k in ["branch", "stack", "mul", "carry", "clmul"]},
                "alias_of": names[target].removeprefix(PROBES) if target != raw else "",
                "external": " ".join(names[e] for e in external),
            }
        )
        listings[probe] = (
            "\n\n".join(f"// {names[r]}\n" + "\n".join(funcs[r].lines) for r in reach)
            + "\n"
        )
    return rows, listings


def table(rows: list[dict], configs: list[str], value) -> list[str]:
    cells = {(r["probe"], r["config"]): value(r) for r in rows}
    probes = sorted({r["probe"] for r in rows}, key=lambda p: (p.split("::")[0], p))
    out = ["| probe | " + " | ".join(configs) + " |", "|---" * (len(configs) + 1) + "|"]
    for p in probes:
        out.append(
            f"| {p} | " + " | ".join(str(cells.get((p, c), "")) for c in configs) + " |"
        )
    return out


def where(rows: list[dict], probe: str, key: str, value: str) -> list[str]:
    """The configurations in which `probe` has `value` for `key`."""
    return sorted(
        {r["config"] for r in rows if r["probe"] == probe and r[key] == value}
    )


def summary(rows: list[dict], configs: list[str]) -> str:
    muls = [r for r in rows if r["probe"].rsplit("::", 1)[-1] in {"mul", "square"}]
    externals = sorted({(r["probe"], r["external"]) for r in rows if r["external"]})
    merged = sorted({(r["probe"], r["alias_of"]) for r in rows if r["alias_of"]})
    lines = [
        "# Assembly probes",
        "",
        "Static counts, not timings: a loop counts once. See asm_report.py.",
        "",
        "## Instructions in each probe's closure (the probe and its crate-local callees)",
        "",
        *table(rows, configs, lambda r: r["closure"]),
        "",
        "## Field multiplies and squares: multiplies / carry ops / carry-less multiplies / stack accesses",
        "",
        *table(
            muls,
            configs,
            lambda r: f"{r['mul']}/{r['carry']}/{r['clmul']}/{r['stack']}",
        ),
        "",
        "## Probes the compiler merged with an identical function",
        "",
    ]
    for probe, target in merged:
        lines.append(
            f"- {probe} = {target} ({', '.join(where(rows, probe, 'alias_of', target))})"
        )
    lines += ["", "## Calls out of the crate", ""]
    for probe, ext in externals:
        lines.append(
            f"- {probe} ({', '.join(where(rows, probe, 'external', ext))}): {ext}"
        )
    return "\n".join(lines) + "\n"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("out", type=Path)
    ap.add_argument("listings", nargs="+", metavar="CONFIG=FILE.s")
    ap.add_argument("--demangler", default="llvm-cxxfilt")
    a = ap.parse_args(argv)
    rows, configs = [], []
    for spec in a.listings:
        config, _, path = spec.partition("=")
        configs.append(config)
        r, listings = probe_rows(
            config, Path(path).read_text(), lambda ns: demangle(ns, a.demangler)
        )
        if not r:
            sys.exit(
                f"asm_report: no probes in {path}; built with --features asm-probes?"
            )
        rows += r
        d = a.out / config
        d.mkdir(parents=True, exist_ok=True)
        for probe, text in listings.items():
            (d / (probe.replace("::", ".") + ".s")).write_text(text)
    keys = ["config", "probe", *COLUMNS, "external"]
    tsv = ["\t".join(keys)] + ["\t".join(str(r[k]) for k in keys) for r in rows]
    (a.out / "probes.tsv").write_text("\n".join(tsv) + "\n")
    (a.out / "summary.md").write_text(summary(rows, configs))


if __name__ == "__main__":
    main()
