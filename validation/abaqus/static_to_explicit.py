"""Abaqus/Explicit quasi-static stand-in for a static deck whose elements Abaqus/Standard does not have (S3RS, S4RS).

    convert: python static_to_explicit.py deck.inp out.inp
    reduce:  python static_to_explicit.py --reduce reference.json deck.inp

convert: every *STATIC step becomes a *DYNAMIC, EXPLICIT step of RAMP + HOLD seconds, its loads (*CLOAD, *DLOAD,
*DSLOAD) follow a smooth step reaching the full value at RAMP and held, and every material gets mass-proportional
damping ALPHA (critical for a fundamental frequency of ALPHA / 2): the model ends at rest in static equilibrium.
reduce: keeps the first and the last frame of the reference, the last one at the end time of the static step, and of
the histories only the initial and final strain and kinetic energies (the increments and iterations are those of
the explicit run).
Runs with the Python of Abaqus (2.7) as well as Python 3.
"""
from __future__ import print_function
import json
import re
import sys

RAMP, HOLD, ALPHA = 0.004, 0.004, 5.6e4


def convert(src, dst):
    lines = open(src).read().splitlines()
    out, i, amp_done = [], 0, False
    while i < len(lines):
        ln = lines[i]
        key = ln.strip().lower()
        if key.startswith('*step') and not amp_done:
            out += ['*Amplitude, name=QSRAMP, definition=SMOOTH STEP',
                    '0., 0., %g, 1., %g, 1.' % (RAMP, RAMP + HOLD)]
            amp_done = True
        if key.startswith('*static'):
            out += ['*Dynamic, Explicit', ', %g' % (RAMP + HOLD)]
            i += 1
            if i < len(lines) and not lines[i].startswith('*'):
                i += 1  # data line of *STATIC
            continue
        if key.startswith('*step'):  # INC is a Standard parameter
            ln = re.sub(r',\s*inc\s*=\s*[^,]*', '', ln, flags=re.I)
        if re.match(r'\*(cload|dload|dsload)\b', key) and 'amplitude' not in key:
            ln += ', amplitude=QSRAMP'
        out.append(ln)
        if key.startswith('*material'):
            out += ['*Damping, alpha=%g' % ALPHA]
        i += 1
    out.insert(1, '** quasi-static Abaqus/Explicit stand-in (static_to_explicit.py): ramp %g s, hold %g s, '
                  'mass-proportional damping %g' % (RAMP, HOLD, ALPHA))
    open(dst, 'w').write('\n'.join(out) + '\n')


def static_end_time(deck):
    """Total time at the end of the static steps (step times from the *STATIC data lines, default 1)."""
    lines, t = open(deck).read().splitlines(), 0.0
    for k, ln in enumerate(lines):
        if ln.strip().lower().startswith('*static'):
            data = lines[k + 1] if k + 1 < len(lines) and not lines[k + 1].startswith('*') else ''
            f = [x.strip() for x in data.split(',')]
            t += float(f[1]) if len(f) > 1 and f[1] else 1.0
    return t or 1.0


def reduce(ref_path, deck):
    d = json.load(open(ref_path))
    fr = d['frames']
    if not fr:
        raise SystemExit('static_to_explicit: no frame in %s (explicit run failed?)' % ref_path)
    last = dict(fr[-1])
    t_end = static_end_time(deck)
    last['t'] = t_end
    d['frames'] = [fr[0], last]
    # energies: initial and final state (strain and kinetic energy; the external work of the explicit run also
    # holds the damping losses)
    d['history'] = dict((k, [[0.0, d['history'][k][0][1]], [t_end, d['history'][k][-1][1]]])
                        for k in ('ALLIE', 'ALLKE') if d['history'].get(k))
    for k in ('rf_hist', 'mon_hist'):
        d.pop(k, None)
    d['abaqus']['note'] = ('Abaqus/Explicit quasi-static (element not available in Abaqus/Standard): ramp %g s, '
                           'hold %g s, mass-proportional damping %g; final state' % (RAMP, HOLD, ALPHA))
    json.dump(d, open(ref_path, 'w'))


if __name__ == '__main__':
    if sys.argv[1] == '--reduce':
        reduce(sys.argv[2], sys.argv[3])
    else:
        convert(sys.argv[1], sys.argv[2])
