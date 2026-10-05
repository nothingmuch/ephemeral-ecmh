import asm_report as ar
import pytest

# Trimmed from rustc's output: ELF aarch64 (Android), Mach-O aarch64 with
# Apple's arrangement suffixes and ; comments, ELF x86-64 in Intel syntax.
ELF_ARM = """\
	.section	.text.mul,"ax",@progbits
	.globl	probe_mul
	.type	probe_mul,@function
probe_mul:
	.cfi_startproc
	str	x30, [sp, #-16]!   // spill
	mul	x8, x0, x1
	umulh	x9, x0, x1
	adds	x8, x8, x9
	cbz	x8, .LBB0_2
	bl	helper
.LBB0_2:
	ldr	x30, [sp], #16
	b	memcpy
.Lfunc_end0:
	.size	probe_mul, .Lfunc_end0-probe_mul
	.cfi_endproc

helper:
	.cfi_startproc
	and	x0, x0, #0x1111
	blr	x8
	ret
	.cfi_endproc
probe_sub = probe_mul
"""

MACHO_ARM = """\
	.globl	__probe_clmul
	.p2align	2
__probe_clmul:
	.cfi_startproc
	ldp	d0, d1, [x1]
	pmull.1q	v4, v2, v0  ; low halves
	pmull2.2d	v5, v3, v1
	eor.16b	v0, v0, v1
	b.ne	LBB1_1
LBB1_1:
	adrp	x9, _gotfn@GOTPAGE
	ldr	x9, [x9, _gotfn@GOTPAGEOFF]
	blr	x9
	ret
	.cfi_endproc
"""

X86 = """\
probe_x86:
	.cfi_startproc
	push	rbx
	mulx	rax, rcx, rdx   # high half
	adc	rax, 0
	vpclmulqdq	xmm0, xmm1, xmm2, 17
	mov	qword ptr [rsp + 8], rax
	jne	.LBB2_1
	call	qword ptr [rip + far@GOTPCREL]
	mov	rbx, qword ptr [rip + far2@GOTPCREL]
	call	rbx
.LBB2_1:
	pop	rbx
	jmp	near@PLT
	.cfi_endproc
"""


def names(ns):
    return {
        n: "ephemeral_ecmh::asm_probes::" + n.removeprefix("__").removeprefix("probe_")
        if "probe" in n
        else n
        for n in ns
    }


def test_parse_functions_and_aliases():
    funcs, aliases = ar.parse(ELF_ARM, "aarch64")
    assert list(funcs) == ["probe_mul", "helper"]
    assert [m for m, _ in funcs["probe_mul"].insns] == [
        "str",
        "mul",
        "umulh",
        "adds",
        "cbz",
        "bl",
        "ldr",
        "b",
    ]
    assert funcs["probe_mul"].insns[0] == ("str", "x30, [sp, #-16]!")
    assert aliases == {"probe_sub": "probe_mul"}


def test_aliases_resolve_through_chains():
    _, aliases = ar.parse(
        "a = b\nb = c\nc:\n\t.cfi_startproc\n\tret\n\t.cfi_endproc\n", "aarch64"
    )
    assert aliases == {"a": "c", "b": "c"}


def test_targets_skip_local_labels_and_registers():
    funcs, _ = ar.parse(ELF_ARM, "aarch64")
    assert ar.targets(funcs["probe_mul"].insns, "aarch64") == ["helper", "memcpy"]
    assert ar.targets(funcs["helper"].insns, "aarch64") == []
    funcs, _ = ar.parse(X86, "x86_64")
    assert ar.targets(funcs["probe_x86"].insns, "x86_64") == ["far", "far2", "near"]
    funcs, _ = ar.parse(MACHO_ARM, "aarch64")
    assert ar.targets(funcs["__probe_clmul"].insns, "aarch64") == ["_gotfn"]
    elf_got = [
        ("adrp", "x8, :got:gotfn"),
        ("ldr", "x8, [x8, :got_lo12:gotfn]"),
        ("blr", "x8"),
    ]
    assert ar.targets(elf_got, "aarch64") == ["gotfn"]


def test_stats_aarch64():
    funcs, _ = ar.parse(ELF_ARM, "aarch64")
    s = ar.stats(funcs["probe_mul"].insns, "aarch64")
    assert s == {"insns": 8, "mul": 2, "carry": 1, "branch": 1, "call": 1, "stack": 2}


def test_stats_apple_suffixes():
    funcs, _ = ar.parse(MACHO_ARM, "aarch64")
    s = ar.stats(funcs["__probe_clmul"].insns, "aarch64")
    assert s["clmul"] == 2 and s["branch"] == 1 and s["insns"] == 9


def test_stats_x86():
    funcs, _ = ar.parse(X86, "x86_64")
    s = ar.stats(funcs["probe_x86"].insns, "x86_64")
    assert s == {
        "insns": 11,
        "mul": 1,
        "carry": 1,
        "clmul": 1,
        "stack": 3,
        "branch": 1,
        "call": 2,
    }


def test_rows_follow_closure_and_aliases():
    rows, listings = ar.probe_rows("aarch64-linux-android/x", ELF_ARM, names)
    by = {r["probe"]: r for r in rows}
    assert by["mul"]["insns"] == 8
    assert by["mul"]["closure"] == 11
    assert by["mul"]["callees"] == 1
    assert by["mul"]["external"] == "memcpy"
    assert by["sub"]["closure"] == 11 and by["sub"]["alias_of"] == "mul"
    assert "helper:" in listings["mul"]
    assert "helper" not in by


def test_arch_must_be_known():
    with pytest.raises(ValueError):
        ar.arch_of("riscv64gc-unknown-linux-gnu/generic")


def test_main_writes_tables(tmp_path, monkeypatch):
    s = tmp_path / "a.s"
    s.write_text(ELF_ARM)
    monkeypatch.setattr(ar, "demangle", lambda ns, _: names(ns))
    ar.main([str(tmp_path / "out"), f"aarch64-linux-android/x={s}"])
    tsv = (tmp_path / "out" / "probes.tsv").read_text().splitlines()
    assert tsv[0].split("\t")[:4] == ["config", "probe", "insns", "closure"]
    assert len(tsv) == 3
    md = (tmp_path / "out" / "summary.md").read_text()
    assert "| mul | 11 |" in md
    assert "| mul | 2/1/0/2 |" in md
    assert "- sub = mul (aarch64-linux-android/x)" in md
    assert (tmp_path / "out" / "aarch64-linux-android/x" / "mul.s").exists()


def test_main_rejects_listings_without_probes(tmp_path, monkeypatch):
    s = tmp_path / "a.s"
    s.write_text("f:\n\t.cfi_startproc\n\tret\n\t.cfi_endproc\n")
    monkeypatch.setattr(ar, "demangle", lambda ns, _: {n: n for n in ns})
    with pytest.raises(SystemExit):
        ar.main([str(tmp_path / "out"), f"x86_64-unknown-linux-gnu/x={s}"])
