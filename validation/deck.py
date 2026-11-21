"""Building blocks of the validation decks: materials, flat input decks, meshes, steps, case registry.

Every case is written in its native procedure; register() also produces its twin in the other procedure when the
generator runs in twin mode (MODE['twin']): a static case becomes a quasi-static Abaqus/Explicit-style step (loads and
boundary conditions follow a smooth step), an explicit case becomes an implicit dynamic step (*DYNAMIC).
Units: mm, N, tonne, s, MPa, mJ.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cases')

MATERIALS = {
    'STEEL': '*Material, name=STEEL\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n*Expansion\n1.2e-05\n',
    'STEEL0': '*Material, name=STEEL0\n*Density\n7.85e-09\n*Elastic\n210000., 0.\n',
    'STEEL_PL': ('*Material, name=STEEL_PL\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n'
                 '*Plastic\n250., 0.\n300., 0.05\n340., 0.2\n360., 1.\n'),
    'ALU': '*Material, name=ALU\n*Density\n2.7e-09\n*Elastic\n70000., 0.33\n',
    'RUBBER': '*Material, name=RUBBER\n*Density\n1.1e-09\n*Hyperelastic, neo hooke\n0.5, 0.04\n',
    'FUNG': ('*Material, name=FUNG\n*Density\n1.06e-09\n'
             '*Anisotropic Hyperelastic, fung-anisotropic, local directions=2\n'
             '8., 0., 100., 0., 0., 8., 0., 0.\n'
             '0., 12., 0., 0., 0., 0., 8., 0.\n'
             '0., 0., 0., 0., 12., 2e-05, 2e-05\n'),
    'STEEL49': '*Material, name=STEEL49\n*Density\n7.85e-09\n*Elastic\n210000., 0.49\n',
    'NH45': '*Material, name=NH45\n*Density\n1.1e-09\n*Hyperelastic, neo hooke\n100., 0.001\n',
    'NH499': '*Material, name=NH499\n*Density\n1.1e-09\n*Hyperelastic, neo hooke\n100., 2e-05\n',
    'HEAVY': '*Material, name=HEAVY\n*Density\n7.85e-06\n*Elastic\n210000., 0.3\n',
    'SOFT7': '*Material, name=SOFT7\n*Density\n7.85e-09\n*Elastic\n7000., 0.33\n',
    'LIGHT7': '*Material, name=LIGHT7\n*Density\n7.85e-10\n*Elastic\n70000., 0.33\n',
    'STEEL_DAMP': '*Material, name=STEEL_DAMP\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n*Damping, alpha=200.\n',
    'STEEL_TDEP': ('*Material, name=STEEL_TDEP\n*Density\n7.85e-09\n*Elastic\n210000., 0.3, 0.\n150000., 0.3, 200.\n'
                   '*Expansion\n1.2e-05\n'),
    'MOONEY': '*Material, name=MOONEY\n*Density\n1.1e-09\n*Hyperelastic, mooney-rivlin\n0.3, 0.2, 0.04\n',
    'STEEL_KIN': ('*Material, name=STEEL_KIN\n*Density\n7.85e-09\n*Elastic\n210000., 0.3\n'
                  '*Plastic, hardening=KINEMATIC\n250., 0.\n350., 0.1\n'),
    'VISCO': ('*Material, name=VISCO\n*Density\n1.1e-09\n*Elastic\n10., 0.45\n'
              '*Viscoelastic, time=PRONY\n0.5, 0., 0.001\n'),
}

TOL = {  # relative to the peak reference value (U also for rotations, RF also for moments, S for element stresses)
    'exact': {'U': 0.001, 'RF': 0.001, 'S': 0.001, 'E': 0.001},    # closed form = exact solution of the discrete model
    'theory': {'U': 0.001, 'RF': 0.001, 'S': 0.002, 'E': 0.001},   # closed form, discretisation error < 0.05 %
    'linear': {'U': 0.02, 'RF': 0.02, 'S': 0.05, 'E': 0.03},       # Abaqus references
    'nonlinear': {'U': 0.04, 'RF': 0.04, 'S': 0.08, 'E': 0.05},
    'dynamic': {'U': 0.06, 'RF': 0.10, 'S': 0.12, 'E': 0.06},
    'contact': {'U': 0.15, 'RF': 0.20, 'S': 0.25, 'E': 0.15},
    'published': {'U': 0.02, 'RF': 0.02, 'S': 0.03, 'E': 0.03},    # published values: tolerance of each check
}

MODE = {'twin': False}
CURRENT = {'deck': None}
QS_DURATION = {}     # case name -> duration of its quasi-static explicit twin (s)
CASES = []
KNOWN = {}           # case name -> documented limitation (KNOWN FAIL)
KNOWN_QUANTITIES = {}  # case name -> (quantities, note): documented difference on some quantities only


def fmt(v):
    return ('%.10g' % v) if abs(v) > 1e-14 else '0.'


def chunks(seq, n=16):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def near(a, b):
    return abs(a - b) < 1e-9


class Deck:
    """Flat input deck (no parts), nodes and elements numbered from 1."""

    def __init__(self, name, title):
        self.name, self.title = name, title
        self.nodes = {}            # id -> (x, y, z)
        self.blocks = []           # (type, elset, [(id, conn)])
        self.nsets, self.elsets = {}, {}
        self.model = []            # model-level text (sections, orientations, surfaces, interactions, constraints)
        self.materials = set()
        self.initial = []          # initial conditions text
        self.steps = []
        self.next_node, self.next_elem = 1, 1
        self.bc_state = {}         # last prescribed value per (target, dof), for quasi-static twins
        self.namp = 0
        CURRENT['deck'] = self

    def node(self, x, y, z=0.0):
        i = self.next_node
        self.nodes[i] = (x, y, z)
        self.next_node += 1
        return i

    def elements(self, etype, elset, conns):
        lst = []
        for c in conns:
            lst.append((self.next_elem, c))
            self.next_elem += 1
        self.blocks.append((etype, elset, lst))
        self.elsets.setdefault(elset, []).extend(e for e, _ in lst)
        return [e for e, _ in lst]

    def nset(self, name, ids):
        self.nsets[name] = sorted(set(ids))

    def select(self, pred):
        return [i for i, p in self.nodes.items() if pred(*p)]

    def text(self):
        out = ['*Heading', self.title, '*Node']
        out += ['%d, %s, %s, %s' % (i, fmt(p[0]), fmt(p[1]), fmt(p[2])) for i, p in sorted(self.nodes.items())]
        for etype, elset, lst in self.blocks:
            out.append('*Element, type=%s, elset=%s' % (etype, elset))
            out += ['%d, %s' % (e, ', '.join(str(n) for n in c)) for e, c in lst]
        for name, ids in self.nsets.items():
            out.append('*Nset, nset=%s' % name)
            out += [', '.join(str(i) for i in ch) for ch in chunks(ids)]
        for name, ids in self.elsets.items():
            out.append('*Elset, elset=%s' % name)
            out += [', '.join(str(i) for i in ch) for ch in chunks(ids)]
        out += self.model
        for m in sorted(self.materials):
            out.append(MATERIALS[m].rstrip('\n'))
        out += self.initial
        out += self.steps
        return '\n'.join(out) + '\n'


# ------------------------------------------------------------------------------------------------ meshes
def hex_block(d, nx, ny, nz, lx, ly, lz, x0=(0.0, 0.0, 0.0)):
    ids = {}
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ids[i, j, k] = d.node(x0[0] + lx * i / nx, x0[1] + ly * j / ny, x0[2] + lz * k / nz)
    hexes = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                hexes.append([ids[i, j, k], ids[i + 1, j, k], ids[i + 1, j + 1, k], ids[i, j + 1, k],
                              ids[i, j, k + 1], ids[i + 1, j, k + 1], ids[i + 1, j + 1, k + 1], ids[i, j + 1, k + 1]])
    return hexes


def vol6(d, c):
    a, b, e, f = (d.nodes[n] for n in c[:4])
    u = [b[i] - a[i] for i in range(3)]
    v = [e[i] - a[i] for i in range(3)]
    w = [f[i] - a[i] for i in range(3)]
    return (u[1] * v[2] - u[2] * v[1]) * w[0] + (u[2] * v[0] - u[0] * v[2]) * w[1] + (u[0] * v[1] - u[1] * v[0]) * w[2]


def tets_from_hexes(d, hexes):
    """Conforming split of each hexahedron into 6 tetrahedra around the diagonal 1-7 (positive volumes)."""
    ring = [(1, 2), (2, 3), (3, 7), (7, 4), (4, 5), (5, 1)]  # vertices around the diagonal 0-6
    tets = []
    for h in hexes:
        for a, b in ring:
            t = [h[0], h[a], h[b], h[6]]
            if vol6(d, t) < 0:
                t[1], t[2] = t[2], t[1]
            tets.append(t)
    return tets


def quad_plate(d, nx, ny, lx, ly, x0=(0.0, 0.0, 0.0)):
    ids = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            ids[i, j] = d.node(x0[0] + lx * i / nx, x0[1] + ly * j / ny, x0[2])
    quads = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]] for j in range(ny) for i in range(nx)]
    return quads, ids


def tris_from_quads(quads):
    tris = []
    for k, q in enumerate(quads):
        if k % 2 == 0:
            tris += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
        else:
            tris += [[q[0], q[1], q[3]], [q[1], q[2], q[3]]]
    return tris


def beam_line(d, n, length, x0=(0.0, 0.0, 0.0), axis=0):
    ids = []
    for i in range(n + 1):
        p = list(x0)
        p[axis] += length * i / n
        ids.append(d.node(*p))
    return [[ids[i], ids[i + 1]] for i in range(n)], ids


HEX_FACES = {1: (0, 1, 2, 3), 2: (4, 5, 6, 7), 3: (0, 1, 5, 4), 4: (1, 2, 6, 5), 5: (2, 3, 7, 6), 6: (3, 0, 4, 7)}
TET_FACES = {1: (0, 1, 2), 2: (0, 3, 1), 3: (1, 3, 2), 4: (2, 3, 0)}
WEDGE_FACES = {1: (0, 1, 2), 2: (3, 5, 4), 3: (0, 3, 4, 1), 4: (1, 4, 5, 2), 5: (2, 5, 3, 0)}


def faces_on(d, etype, elist, conns, pred):
    """(element, face) pairs whose corner nodes all satisfy pred."""
    table = TET_FACES if etype.startswith('C3D4') or etype.startswith('C3D10') else (
        WEDGE_FACES if etype.startswith('C3D6') else HEX_FACES)
    res = []
    for e, c in zip(elist, conns):
        for f, loc in table.items():
            if all(pred(*d.nodes[c[i]]) for i in loc):
                res.append((e, f))
    return res


def edge_shares(d, nodes):
    """Consistent nodal shares of a load uniformly distributed along a straight edge of equal elements."""
    nodes = sorted(nodes, key=lambda n: d.nodes[n][1])
    w = [0.5 if i in (0, len(nodes) - 1) else 1.0 for i in range(len(nodes))]
    return [(n, wi / sum(w)) for n, wi in zip(nodes, w)]


def boundary_stretch(d, nodes, e1, e2):
    """Displacements of a homogeneous in-plane stretch (1 + e1, 1 + e2) on the given nodes."""
    out = []
    for n in sorted(set(nodes)):
        x, y, z = d.nodes[n]
        out += ['%d, 1, 1, %s' % (n, fmt(e1 * x)), '%d, 2, 2, %s' % (n, fmt(e2 * y))]
    return out


# ------------------------------------------------------------------------------------------------- steps
STRUCTURAL = ('C3D', 'S4', 'S3', 'B3', 'M3D', 'T3D')
PLASTIC_MATS = ('STEEL_PL', 'STEEL_KIN', 'STEEL_RATE', 'STEEL_DMG', 'COPPER')


def monitored_elements(d):
    """First, middle and last element of the first structural element block (element output set MONE)."""
    for etype, elset, lst in d.blocks:
        if etype.upper().startswith(STRUCTURAL):
            ids = [e for e, _ in lst]
            return etype.upper(), sorted(set([ids[0], ids[len(ids) // 2], ids[-1]]))
    return None, []


def element_output(d):
    etype, mone = monitored_elements(d)
    if not mone:
        return []
    plastic = any(m in d.materials for m in PLASTIC_MATS)
    if etype.startswith('B3'):
        var = 'SF'
    elif etype.startswith(('S4', 'S3')):
        var = 'S, SF'
    else:
        var = 'S'
    if plastic and not etype.startswith('B3'):
        var += ', PEEQ'
    return ['*Element Output, elset=MONE', var]


def output_block(standard, rot, nframes, duration=None, time_marks=None):
    """Field and history output. Field frames at exact times (TIME MARKS=YES) for Abaqus/Standard; Abaqus/Explicit
    by default writes them at the end of the first increment past each output time (time_marks=True to force exact
    times, e.g. for a closed-form reference)."""
    d = CURRENT['deck']
    nv = 'U, RF' + (', UR, RM' if rot else '') + ('' if standard else ', V')
    exact = standard if time_marks is None else time_marks
    s = ['*Output, field, number interval=%d%s' % (nframes, ', time marks=YES' if exact else ''),
         '*Node Output', nv] + element_output(d)
    s += ['*Output, history, variable=PRESELECT' + ('' if standard else ', time interval=%s' % fmt(duration / 200))]
    if not standard:
        s += ['*Output, history, time interval=%s' % fmt(duration / 200), '*Incrementation Output', 'DT, DMASS']
    return s


LOAD_KW = ('*boundary', '*cload', '*dsload', '*dload', '*temperature', '*field')
CONTACT_KW = ('*contact', '*contact inclusions', '*contact property assignment', '*contact pair')


def _smooth_amplitude(d, dur, a0=0.0):
    d.namp += 1
    name = 'QS%d' % d.namp
    d.model += ['*Amplitude, name=%s, definition=SMOOTH STEP' % name, '0., %s, %s, 1.' % (fmt(a0), fmt(dur))]
    return name


def quasi_static_body(d, body, dur):
    """Loads and prescribed values of a static step, made to follow a smooth step over the explicit step. A value
    changed from a non-zero previous value starts from it (amplitude from previous/new to 1)."""
    out, kw, head = [], None, None
    base = None
    for line in body:
        if line.startswith('*'):
            low = line.lower()
            kw = low.split(',')[0].strip()
            head = line
            if kw in LOAD_KW and 'amplitude=' not in low:
                if base is None:
                    base = _smooth_amplitude(d, dur)
                out.append(line + ', amplitude=%s' % base)
            else:
                out.append(line)
            continue
        p = [v.strip() for v in line.split(',')]
        if '<' in line:  # parameter reference: value unknown here, no continuation amplitude
            out.append(line)
            continue
        if kw in ('*boundary', '*cload', '*temperature', '*field') and 'amplitude=' not in head.lower():
            if kw == '*boundary':
                key, val = (p[0], p[1], p[2] if len(p) > 2 else p[1]), (float(p[3]) if len(p) > 3 and p[3] else 0.0)
            elif kw == '*cload':
                key, val = (kw, p[0], p[1]), float(p[2]) if len(p) > 2 and p[2] else 0.0
            else:
                key, val = (kw, p[0]), float(p[-1]) if p[-1] else 0.0
            prev = d.bc_state.get(key, 0.0)
            d.bc_state[key] = val
            if prev != 0.0 and val != 0.0 and abs(prev - val) > 1e-14 * abs(val):
                amp = _smooth_amplitude(d, dur, prev / val)
                out += [head + ', amplitude=%s' % amp, line, head + ', amplitude=%s' % base]
                continue
        out.append(line)
    return out


def static_step(name, body, nlgeom, rot=False, nframes=10, inc=0.1, maxinc=None):
    d = CURRENT['deck']
    if MODE['twin']:  # quasi-static explicit twin
        dur = QS_DURATION.get(d.name, 0.01)
        head = ['*Step, name=%s, nlgeom=YES' % name, '*Dynamic, Explicit', ', %s' % fmt(dur)]
        return head + quasi_static_body(d, body, dur) + output_block(False, rot, nframes, dur) + ['*End Step']
    if not nlgeom and nframes == 1:
        head = ['*Step, name=%s, nlgeom=NO' % name, '*Static', '1., 1.']
    else:
        mx = maxinc if maxinc else 1.0 / nframes
        head = ['*Step, name=%s, nlgeom=%s, inc=10000' % (name, 'YES' if nlgeom else 'NO'), '*Static',
                '%s, 1., 1e-08, %s' % (fmt(min(inc, mx)), fmt(mx))]
    return head + body + output_block(True, rot, nframes) + ['*End Step']


def explicit_step(name, duration, body, rot=False, nframes=20, extra=None, fixed_dt=None, time_marks=None):
    d = CURRENT['deck']
    extra = extra or []
    if MODE['twin']:  # implicit dynamic twin; contact definitions move to the model data (Abaqus/Standard)
        kept, kw = [], None
        for line in body:
            if line.startswith('*'):
                kw = line.lower().split(',')[0].strip()
            (d.model if kw in CONTACT_KW else kept).append(line)
        head = ['*Step, name=%s, nlgeom=YES, inc=1000000' % name, '*Dynamic, application=MODERATE DISSIPATION',
                '%s, %s, %s, %s' % (fmt(duration / 1000), fmt(duration), fmt(duration * 1e-9), fmt(duration / 100))]
        return head + kept + output_block(True, rot, nframes) + ['*End Step']
    head = ['*Step, name=%s, nlgeom=YES' % name] + (
        ['*Dynamic, Explicit, direct user control', '%s, %s' % (fmt(fixed_dt), fmt(duration))] if fixed_dt else
        ['*Dynamic, Explicit', ', %s' % fmt(duration)]) + extra
    return head + body + output_block(False, rot, nframes, duration, time_marks) + ['*End Step']


# ---------------------------------------------------------------------------------------------- registry
def register(d, procedure, element, mode, description, mon, fix, rot, tol, analysis='nonlinear', known=None,
             reference='abaqus', supported=True, reason=None, twin=True, category='element', group=None,
             literature=None, energies=None):
    """Writes the deck and its manifest entry. In twin mode, writes the twin instead (Abaqus-referenced cases only)."""
    d.nset('MON', mon)
    if fix:
        d.nset('FIX', fix)
    etype, mone = monitored_elements(d)
    if mone:
        d.elsets['MONE'] = mone
    name = d.name
    if MODE['twin']:
        if reference != 'abaqus' or not twin:
            return
        if procedure == 'static':
            name, procedure, analysis = name + '_xpl', 'explicit', 'quasi-static explicit'
            tol = 'contact' if tol == 'contact' else 'dynamic'
        else:
            name, procedure, analysis = name + '_std', 'implicit dynamic', 'implicit dynamic'
            supported, reason = False, 'implicit dynamics (*DYNAMIC) not implemented'
        known = None
    kq = KNOWN_QUANTITIES.get(name, ([], None))
    entry = {'name': name, 'category': category, 'group': group or element, 'procedure': procedure,
             'element': element, 'mode': mode, 'description': description, 'mon': sorted(set(mon)),
             'fix': sorted(set(fix)), 'mone': mone, 'mone_type': etype, 'rot': rot, 'tol': TOL[tol],
             'analysis': analysis,
             'known': known or KNOWN.get(name), 'reference': reference, 'known_quantities': kq[0],
             'known_quantities_note': kq[1], 'supported': supported, 'unsupported_reason': reason,
             'twin_of': d.name if MODE['twin'] else None, 'literature': literature, 'energies': energies,
             'deferred': (group == 'application') and ('application validation case, addressed once the simple '
                                                       'tests are validated')}
    CASES.append(entry)
    text = d.text()
    if getattr(d, 'import_only', False):  # model imported from a previous job: no mesh, no sets
        text = '*Heading\n%s\n%s\n' % (d.title, '\n'.join(d.steps))
    elif getattr(d, 'include_mesh', False):  # mesh and sets in an included file, values from *PARAMETER
        lines = text.split('\n')
        k0 = lines.index('*Node')
        k1 = next(k for k, l in enumerate(lines) if l.startswith('*Solid Section') or l.startswith('*Shell Section'))
        mesh = lines[k0:k1] + ['*Elset, elset=EGEN, generate', '1, %d, 1' % (d.next_elem - 1)]
        with open(os.path.join(OUT, name + '_mesh.inc'), 'w') as f:
            f.write('\n'.join(mesh) + '\n')
        text = '\n'.join(lines[:k0] + ['*Parameter', 'YOUNG = 210000.', 'NU = 0.3', 'TIPU = 0.04',
                                       '*Include, input=%s_mesh.inc' % name] + lines[k1:])
    if procedure == 'explicit' and reference == 'abaqus' and fix:
        text = rf_history(text, rot)
    with open(os.path.join(OUT, name + '.inp'), 'w') as f:
        f.write(text.replace('*Heading\n%s' % d.title, '*Heading\n%s%s' % (
            d.title, '' if not MODE['twin'] else (' [quasi-static explicit twin]' if name.endswith('_xpl')
                                                   else ' [implicit dynamic twin]'))))


def rf_history(text, rot=False):
    """Reaction forces of the FIX nodes and motion of the MON nodes at every increment of each explicit step
    (history output): reactions and rotations of dynamic cases are compared after a moving average over a few tens of
    increments (see compare.py)."""
    out, explicit = [], False
    for line in text.split('\n'):
        low = line.lower()
        if low.startswith('*step'):
            explicit = False
        elif low.startswith('*dynamic') and 'explicit' in low:
            explicit = True
        elif low.startswith('*end step') and explicit:
            out += ['*Output, history, frequency=1', '*Node Output, nset=FIX', 'RF1, RF2, RF3',
                    '*Node Output, nset=MON', 'U1, U2, U3' + (', UR1, UR2, UR3' if rot else '')]
        out.append(line)
    return '\n'.join(out)


def write_manifest():
    with open(os.path.join(OUT, 'manifest.json'), 'w') as f:
        json.dump(CASES, f, indent=1)
