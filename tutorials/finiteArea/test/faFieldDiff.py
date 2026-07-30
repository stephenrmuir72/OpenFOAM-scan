#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# faFieldDiff.py  --  finiteArea tutorial parity/regression comparator
#
# Compares the internalField of two OpenFOAM fields (scalar or vector) and
# reports Linf / L2 / relative-L2 norms of the difference, gating on explicit
# tolerances.  Used by the finiteArea tutorial CI harness (tutorials/
# finiteArea/Alltest) to enforce parity of a migrated finiteArea case against
# a stored reference (see that harness for how references are generated).
#
# The comparison is index-wise: it assumes both fields are defined on the same
# area mesh with the same face ordering.  This holds when both meshes are built
# from the same faMeshDefinition/parent mesh (legacy makeFaMesh and makeFaMesh
# both derive faceLabels by iterating the selected polyPatch faces in the same
# order).  As a guard the comparator refuses to compare fields of unequal
# length.
#
# Usage:
#   faFieldDiff.py REF NEW [--abs-tol A] [--rel-tol R] [--name TAG]
#
# Exit status 0 => within tolerance (PASS), 1 => exceeds tolerance (FAIL),
# 2 => usage/parse error.
# -----------------------------------------------------------------------------

import re
import sys
import math
import argparse


def _read_tokens(path):
    with open(path, "r") as f:
        return f.read()


def parse_internal_field(path):
    """Return (kind, data) where kind is 'scalar' or 'vector' and data is a
    list of floats (scalar) or list of 3-tuples (vector)."""
    txt = _read_tokens(path)

    # Locate the internalField entry.
    m = re.search(r"internalField\s+(.*?);", txt, re.DOTALL)
    if not m:
        raise ValueError(f"{path}: no internalField entry found")
    body = m.group(1).strip()

    is_vector = "List<vector>" in body or body.startswith("uniform (")

    if body.startswith("uniform"):
        # uniform <value> : need the mesh size to expand -> not supported here,
        # callers compare nonuniform solution fields.  Return single value.
        val = body[len("uniform"):].strip()
        if val.startswith("("):
            comps = [float(x) for x in val.strip("()").split()]
            return ("vector", [tuple(comps)])
        return ("scalar", [float(val)])

    # nonuniform List<...> N ( ... )
    lp = body.find("(")
    rp = body.rfind(")")
    if lp < 0 or rp < 0:
        raise ValueError(f"{path}: malformed nonuniform list")
    inner = body[lp + 1:rp]

    if is_vector:
        # entries like "(a b c)"
        vecs = re.findall(r"\(([^()]*)\)", inner)
        data = []
        for v in vecs:
            comps = [float(x) for x in v.split()]
            data.append(tuple(comps))
        return ("vector", data)
    else:
        data = [float(x) for x in inner.split()]
        return ("scalar", data)


def norms(ref, new, kind):
    n = len(ref)
    linf = 0.0
    l2 = 0.0
    refl2 = 0.0
    if kind == "scalar":
        for a, b in zip(ref, new):
            d = abs(b - a)
            linf = max(linf, d)
            l2 += d * d
            refl2 += a * a
    else:
        for a, b in zip(ref, new):
            d2 = sum((b[i] - a[i]) ** 2 for i in range(len(a)))
            d = math.sqrt(d2)
            linf = max(linf, d)
            l2 += d2
            refl2 += sum(a[i] * a[i] for i in range(len(a)))
    l2 = math.sqrt(l2 / max(n, 1))
    refl2 = math.sqrt(refl2 / max(n, 1))
    rel = l2 / refl2 if refl2 > 0 else (0.0 if l2 == 0 else float("inf"))
    return linf, l2, rel


def main():
    ap = argparse.ArgumentParser(description="finiteArea field parity diff")
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--abs-tol", type=float, default=1e-8)
    ap.add_argument("--rel-tol", type=float, default=1e-6)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()

    tag = args.name or f"{args.ref} vs {args.new}"

    try:
        kref, dref = parse_internal_field(args.ref)
        knew, dnew = parse_internal_field(args.new)
    except Exception as e:  # noqa: BLE001
        print(f"[faFieldDiff] PARSE-ERROR {tag}: {e}", file=sys.stderr)
        return 2

    if kref != knew:
        print(f"[faFieldDiff] TYPE-MISMATCH {tag}: {kref} vs {knew}",
              file=sys.stderr)
        return 2
    if len(dref) != len(dnew):
        print(f"[faFieldDiff] SIZE-MISMATCH {tag}: "
              f"{len(dref)} vs {len(dnew)} (different mesh ordering/size)",
              file=sys.stderr)
        return 2

    linf, l2, rel = norms(dref, dnew, kref)
    ok = (linf <= args.abs_tol) or (rel <= args.rel_tol)
    status = "PASS" if ok else "FAIL"
    print(f"[faFieldDiff] {status} {tag}: n={len(dref)} kind={kref} "
          f"Linf={linf:.3e} L2={l2:.3e} relL2={rel:.3e} "
          f"(absTol={args.abs_tol:.1e} relTol={args.rel_tol:.1e})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
