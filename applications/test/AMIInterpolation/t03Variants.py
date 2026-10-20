"""Write T03 cyclicAMI variant cases (blockMeshDict and system/t03Variant).

The neighbour block is placed by M(p) = R(p - c) + c + s, where R is the
rotation (rotational) or identity (translational), c the rotation centre
and s the separation vector. The far displacement is d = R*(0,-1,0).
With the default fixture (90 deg about z, centre 0, s = 0) the generated
vertices reproduce cases/topology/system/blockMeshDict exactly.
"""

from __future__ import annotations

import math
import os

OWNER_VERTS = [
    (1, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0),
    (1, 0, 1), (2, 0, 1), (2, 1, 1), (1, 1, 1),
]


def unit(v):
    m = math.sqrt(sum(x * x for x in v))
    return tuple(x / m for x in v)


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def rotate(v, axis, angle_deg):
    k = unit(axis)
    t = math.radians(angle_deg)
    c, s = math.cos(t), math.sin(t)
    kv = cross(k, v)
    kd = dot(k, v)
    return tuple(v[i] * c + kv[i] * s + k[i] * kd * (1 - c) for i in range(3))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def fmt(v):
    out = []
    for x in v:
        x = round(float(x), 12) + 0.0
        out.append(str(int(x)) if x.is_integer() else repr(x))
    return "(%s)" % " ".join(out)


def write_case(case, transform, axis, centre, angle, sep):
    if transform == "rotational":
        def place(p):
            return add(rotate(sub(p, centre), axis, angle), add(centre, sep))
        far = rotate((0, -1, 0), axis, angle)
    elif transform == "translational":
        def place(p):
            return add(p, sep)
        far = (0, -1, 0)
    else:
        raise ValueError(f"unknown transform: {transform}")

    m = [place(p) for p in OWNER_VERTS]
    nbr = [
        m[0], add(m[0], far), add(m[1], far), m[1],
        m[4], add(m[4], far), add(m[5], far), m[5],
    ]
    verts = OWNER_VERTS + nbr

    if transform == "rotational":
        patch_kw = (
            "transform rotational;\n"
            "        rotationAxis %s;\n"
            "        rotationCentre %s;\n"
            "        rotationAngle %s;\n"
        ) % (fmt(axis), fmt(centre), repr(float(angle)))
    else:
        patch_kw = (
            "transform translational;\n"
            "        separationVector %s;\n"
        ) % fmt(sep)

    block_mesh_dict = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object blockMeshDict;
}

vertices
(
%s
);

blocks
(
    hex (0 1 2 3 4 5 6 7) left (1 1 1) simpleGrading (1 1 1)
    hex (8 9 10 11 12 13 14 15) right (1 1 1) simpleGrading (1 1 1)
);

edges ();

boundary
(
    AMI1
    {
        type cyclicAMI;
        neighbourPatch AMI2;
        %s        lowWeightCorrection 0.5;
        faces ((0 1 5 4));
    }
    AMI2
    {
        type cyclicAMI;
        neighbourPatch AMI1;
        %s        lowWeightCorrection 0.5;
        faces ((8 12 15 11));
    }
);

defaultPatch
{
    name walls;
    type wall;
}

mergePatchPairs ();
""" % (
        "\n".join("    " + fmt(v) for v in verts),
        patch_kw,
        patch_kw,
    )

    variant = """FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object t03Variant;
}

transform %s;
rotationAxis %s;
rotationCentre %s;
rotationAngle %s;
separationVector %s;
""" % (
        transform,
        fmt(axis),
        fmt(centre),
        repr(float(angle)),
        fmt(sep),
    )

    os.makedirs(os.path.join(case, "system"), exist_ok=True)
    with open(os.path.join(case, "system", "blockMeshDict"), "w") as f:
        f.write(block_mesh_dict)
    with open(os.path.join(case, "system", "t03Variant"), "w") as f:
        f.write(variant)

    # Boundary keywords such as rotationAngle are written with
    # writePrecision; full precision keeps the transform exact.
    control_dict = os.path.join(case, "system", "controlDict")
    with open(control_dict, encoding="utf-8") as f:
        lines = [line for line in f.read().splitlines()
                 if "writePrecision" not in line]
    lines.append("writePrecision 17;")
    with open(control_dict, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
