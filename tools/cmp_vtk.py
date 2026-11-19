#!/usr/bin/env python3
"""Compare an asclepios frame (.vtu) with a reference surface exported from Abaqus (legacy ASCII .vtk, POLYDATA).

Correspondence: the reference initial state (legacy .vtk of the undeformed surface) is matched to the asclepios
undeformed coordinates (X - U of the frame) by nearest neighbour; the deformed positions are then compared
node by node. A correspondence-free surface distance (nearest reference point) is also reported.

usage: cmp_vtk.py frame.vtu reference_deformed.vtk reference_initial.vtk
"""
import base64
import re
import struct
import sys

import numpy as np
from scipy.spatial import cKDTree


def read_legacy(path):
    t = open(path, 'rb').read().decode('latin1').split()
    i = t.index('POINTS')
    n = int(t[i + 1])
    return np.array(t[i + 3:i + 3 + 3 * n], float).reshape(n, 3)


def read_vtu(path):
    s = open(path).read()

    def arr(name, typ, nc):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*format="binary"[^>]*>\s*([^<]*)<' % name, s)
        b = base64.b64decode(m.group(1).strip())
        nb = struct.unpack('<I', b[:4])[0]
        return np.frombuffer(b[4:4 + nb], dtype=typ).reshape(-1, nc).astype(float)

    return arr('Points', np.float32, 3), arr('U', np.float32, 3)


def stats(d):
    return 'mean %.3f  median %.3f  p95 %.3f  max %.3f' % (d.mean(), np.median(d), np.percentile(d, 95), d.max())


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)
    x, u = read_vtu(sys.argv[1])
    ref = read_legacy(sys.argv[2])
    ref0 = read_legacy(sys.argv[3])
    d0, idx = cKDTree(x - u).query(ref0)
    if d0.max() > 1e-2:
        print('warning: initial correspondence up to %.3g mm' % d0.max())
    mine = x[idx]
    disp_ref = np.linalg.norm(ref - ref0, axis=1)
    disp_mine = np.linalg.norm(mine - ref0, axis=1)
    err = np.linalg.norm(mine - ref, axis=1)
    surf = cKDTree(ref).query(mine)[0]
    print('nodes compared          %d' % len(ref))
    print('reference displacement  ' + stats(disp_ref))
    print('asclepios displacement    ' + stats(disp_mine))
    print('node-to-node distance   ' + stats(err))
    print('surface distance        ' + stats(surf))
    cosang = np.sum((mine - ref0) * (ref - ref0), axis=1) / np.maximum(disp_mine * disp_ref, 1e-12)
    big = disp_ref > 1.0
    print('direction agreement     mean cos %.3f (nodes moving > 1 mm in the reference)' % cosang[big].mean())


if __name__ == '__main__':
    main()
