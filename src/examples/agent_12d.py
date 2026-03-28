"""
Self-evolving agent to find form field(s) F (no dilaton, at most two forms)
that balance R_{MN} = T_{MN} for:

  ds^2_{12} = H^{-3/4} ds^2_{1,1} + H^{1/2} dz_1^2 + H^{-1/2} dz_2^2 + H^{1/4} ds^2_8

The agent maintains a log file (12d_search.log) where it records:
  - What it tried
  - What it observed
  - What it learned
  - What to try next

Run repeatedly or let it loop -- it picks up from where it left off.
"""

import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import sympy as sp
from itertools import combinations, product
from sugra import (HarmonicFunction, warped_product,
                   FormField, exterior_derivative, hodge_star,
                   form_stress_energy, Verifier)

LOG_FILE = os.path.join(os.path.dirname(__file__), '12d_search.log')

# ============================================================================
# Fixed setup
# ============================================================================

d_wv = 2
d_harm = 8

wv_coords = list(sp.symbols('t x1', real=True))
z1, z2 = sp.symbols('z1 z2', real=True)
harmonic_coords = list(sp.symbols(f'y0:{d_harm}', real=True))

coords = wv_coords + [z1, z2] + harmonic_coords
D = len(coords)

hf = HarmonicFunction(transverse_coords=harmonic_coords)
H_func = sp.Function('H')(hf.r_expr)

metric = warped_product(
    warp_factors     = [H_func**sp.Rational(-3, 4),
                        H_func**sp.Rational(1, 2),
                        H_func**sp.Rational(-1, 2),
                        H_func**sp.Rational(1, 4)],
    block_dims       = [d_wv, 1, 1, d_harm],
    block_signatures = ['lorentzian', 'euclidean', 'euclidean', 'euclidean'],
    coordinates      = coords,
)

# One representative per block (t=wv, z1=torus1, z2=torus2, y0=transverse)
REPS = [(0, 't'), (2, 'z1'), (3, 'z2'), (4, 'y0')]
REP_NAMES = [n for _, n in REPS]

# ============================================================================
# Logging
# ============================================================================

def load_log():
    if os.path.exists(LOG_FILE):
        with open(LOG_FILE, 'r') as f:
            return json.load(f)
    return {
        'R': None,
        'tried': [],
        'insights': [],
        'queue': [],
        'solution': None,
    }

def save_log(log):
    with open(LOG_FILE, 'w') as f:
        json.dump(log, f, indent=2, default=str)

def log_msg(log, msg):
    timestamp = time.strftime('%H:%M:%S')
    entry = f"[{timestamp}] {msg}"
    log.setdefault('messages', []).append(entry)
    print(entry)

# ============================================================================
# Physics helpers
# ============================================================================

def make_F(legs, power, form_type='electric'):
    C = FormField(rank=len(legs), dim=D)
    C[legs] = H_func**power
    F_E = exterior_derivative(C, coords)
    if form_type == 'electric':
        return F_E
    return hodge_star(F_E, metric)

def compute_T_reps(F):
    T_full = form_stress_energy(F, metric)
    T = {}
    for idx, name in REPS:
        T[name] = hf.substitute(sp.cancel(T_full[idx, idx]))
    return T

def extract_coeff_and_power(expr):
    """Given expr of the form c * H'^a * H^b / r^k, extract (c, b).
    Returns (coeff, H_power) or None if not in this form."""
    H, Hp, r = hf.H, hf.Hp, hf.r
    # All our R and T expressions are of the form: coeff * H'^2 * H^n
    # possibly with y0, r factors for the transverse component.
    # Factor out H'^2
    expr = sp.cancel(expr)
    if expr == 0:
        return (0, None)
    # Divide by H'^2
    ratio = sp.cancel(expr / Hp**2)
    # Extract H power
    powers = ratio.as_powers_dict()
    h_power = powers.get(H, 0)
    coeff = sp.cancel(ratio / H**h_power)
    return (coeff, h_power)

# ============================================================================
# Experiment types
# ============================================================================

def experiment_single(legs, power, form_type):
    """Try a single form field, return analysis dict."""
    F = make_F(tuple(legs), power, form_type)
    T = compute_T_reps(F)
    analysis = {}
    for name in REP_NAMES:
        t_decomp = extract_coeff_and_power(T[name])
        r_decomp = extract_coeff_and_power(R_computed[name])
        ratio = sp.cancel(R_computed[name] / T[name]) if T[name] != 0 else 'T=0'
        analysis[name] = {
            'T': str(T[name]),
            'T_coeff': str(t_decomp[0]),
            'T_Hpower': str(t_decomp[1]),
            'R_coeff': str(r_decomp[0]),
            'R_Hpower': str(r_decomp[1]),
            'ratio': str(ratio),
            'H_power_match': str(t_decomp[1]) == str(r_decomp[1]),
        }
    return analysis

def experiment_pair(legs1, p1, type1, legs2, p2, type2):
    """Try sum of two form fields with scalar coefficients.

    Computes T1 and T2 separately, then solves c1^2*T1 + c2^2*T2 = R.
    """
    F1 = make_F(tuple(legs1), p1, type1)
    F2 = make_F(tuple(legs2), p2, type2)
    T1 = compute_T_reps(F1)
    T2 = compute_T_reps(F2)

    c1sq, c2sq = sp.symbols('c1sq c2sq')
    equations = []
    for name in REP_NAMES:
        eq = sp.cancel(c1sq * T1[name] + c2sq * T2[name] - R_computed[name])
        if eq != 0:
            equations.append(eq)

    if not equations:
        return {'status': 'trivial', 'c1sq': 1, 'c2sq': 1}

    try:
        sols = sp.solve(equations, [c1sq, c2sq], dict=True)
    except Exception as e:
        return {'status': 'solve_failed', 'error': str(e)}

    results = []
    for s in sols:
        if c1sq not in s or c2sq not in s:
            continue
        v1, v2 = s[c1sq], s[c2sq]
        # Check: must be constant (no H, r, y dependence)
        free = v1.free_symbols | v2.free_symbols
        bad_syms = {hf.H, hf.Hp, hf.r} | set(hf._y)
        if free & bad_syms:
            results.append({'status': 'H_dependent', 'c1sq': str(v1), 'c2sq': str(v2)})
            continue
        # Check positivity
        if v1.is_negative or v2.is_negative:
            results.append({'status': 'negative', 'c1sq': str(v1), 'c2sq': str(v2)})
            continue
        results.append({'status': 'candidate', 'c1sq': str(v1), 'c2sq': str(v2)})

    return results if results else [{'status': 'no_solution'}]

def verify_pair(legs1, p1, type1, c1, legs2, p2, type2, c2):
    """Full verification of a pair solution."""
    F1 = make_F(tuple(legs1), p1, type1) * c1
    F2 = make_F(tuple(legs2), p2, type2) * c2
    F = F1 + F2
    soln = dict(metric=metric, F=F, Phi=0, alpha=0, coords=coords, hf=hf)
    return Verifier(soln).check()

# ============================================================================
# Queue generators
# ============================================================================

ALL_LEGS = []
for r in range(1, 5):
    for subset in combinations([0, 1, 2, 3], r):
        if 0 in subset:  # must include time
            ALL_LEGS.append(list(subset))

POWERS = [-1, -2, 1, sp.Rational(-1, 2), sp.Rational(1, 2)]
TYPES = ['electric', 'magnetic']

def gen_single_experiments():
    """Generate single-form experiments."""
    for legs in ALL_LEGS:
        for p in POWERS:
            for ftype in TYPES:
                yield {
                    'type': 'single',
                    'legs': legs,
                    'power': str(p),
                    'form_type': ftype,
                }

def gen_pair_experiments(log):
    """Generate pair experiments based on insights from single-form results.

    Strategy: pair forms whose H-power structures complement each other.
    Prioritize pairs where one form matches some blocks and the other matches
    the remaining blocks.
    """
    # Find single-form results with partial H-power matches
    good_singles = []
    for entry in log['tried']:
        if entry.get('exp', {}).get('type') != 'single':
            continue
        analysis = entry.get('result', {})
        if not isinstance(analysis, dict):
            continue
        matched = sum(1 for name in REP_NAMES
                      if isinstance(analysis.get(name), dict)
                      and analysis[name].get('H_power_match'))
        if matched > 0:
            good_singles.append((entry['exp'], matched, analysis))

    # Sort by number of matched components (descending)
    good_singles.sort(key=lambda x: -x[1])

    # Pair them: prefer combining forms that match different components
    seen = set()
    for i, (exp1, m1, a1) in enumerate(good_singles):
        for exp2, m2, a2 in good_singles[i+1:]:
            k = (str(exp1['legs']), exp1['form_type'],
                 str(exp2['legs']), exp2['form_type'])
            if k in seen:
                continue
            seen.add(k)
            # Check if they match complementary components
            matched1 = {n for n in REP_NAMES
                        if isinstance(a1.get(n), dict)
                        and a1[n].get('H_power_match')}
            matched2 = {n for n in REP_NAMES
                        if isinstance(a2.get(n), dict)
                        and a2[n].get('H_power_match')}
            # Prioritize if they cover all components together
            coverage = len(matched1 | matched2)
            yield {
                'type': 'pair',
                'legs1': exp1['legs'],
                'power1': exp1['power'],
                'type1': exp1['form_type'],
                'legs2': exp2['legs'],
                'power2': exp2['power'],
                'type2': exp2['form_type'],
                'coverage': coverage,
            }

    # Also try all pairs with p=-1 if no insights yet
    if not good_singles:
        for i, legs1 in enumerate(ALL_LEGS):
            for legs2 in ALL_LEGS[i:]:
                for t1, t2 in product(TYPES, TYPES):
                    yield {
                        'type': 'pair',
                        'legs1': legs1, 'power1': '-1', 'type1': t1,
                        'legs2': legs2, 'power2': '-1', 'type2': t2,
                    }

# ============================================================================
# Main loop
# ============================================================================

def already_tried(log, exp):
    """Check if this experiment was already run."""
    for entry in log['tried']:
        if entry.get('exp') == exp:
            return True
    return False

def run_agent():
    log = load_log()

    # Phase 0: Compute R if not cached
    global R_computed
    log_msg(log, "Computing Ricci tensor ...")
    R_full = metric.ricci_tensor(simplify_func=sp.cancel)
    R_computed = {}
    for idx, name in REPS:
        R_computed[name] = hf.substitute(sp.cancel(R_full[idx, idx]))
    log['R'] = {k: str(v) for k, v in R_computed.items()}
    log_msg(log, f"R = {log['R']}")
    save_log(log)

    if log.get('solution'):
        log_msg(log, f"Solution already found: {log['solution']}")
        return

    # Phase 1: Single forms
    log_msg(log, "--- Phase 1: Single forms ---")
    for exp in gen_single_experiments():
        if already_tried(log, exp):
            continue

        legs = exp['legs']
        p = sp.S(exp['power'])
        ftype = exp['form_type']
        leg_names = ','.join(str(coords[i]) for i in legs)
        prefix = '' if ftype == 'electric' else '*'
        label = f"{prefix}dC({leg_names}), p={p}"

        log_msg(log, f"Trying: {label}")
        try:
            result = experiment_single(legs, p, ftype)
        except Exception as e:
            result = {'error': str(e)}

        # Learn: check if any components have matching H-powers
        matched = []
        mismatched = []
        for name in REP_NAMES:
            if isinstance(result.get(name), dict):
                if result[name].get('H_power_match'):
                    matched.append(name)
                else:
                    mismatched.append(name)

        insight = f"{label}: H-power matches {matched}, mismatches {mismatched}"
        log_msg(log, f"  -> {insight}")

        log['tried'].append({'exp': exp, 'result': result})
        if insight not in log['insights']:
            log['insights'].append(insight)
        save_log(log)

    # Phase 2: Pairs
    log_msg(log, "--- Phase 2: Pairs of forms ---")
    for exp in gen_pair_experiments(log):
        if already_tried(log, exp):
            continue

        legs1 = exp['legs1']
        p1 = sp.S(exp['power1'])
        t1 = exp['type1']
        legs2 = exp['legs2']
        p2 = sp.S(exp['power2'])
        t2 = exp['type2']

        l1 = ','.join(str(coords[i]) for i in legs1)
        l2 = ','.join(str(coords[i]) for i in legs2)
        pf1 = '' if t1 == 'electric' else '*'
        pf2 = '' if t2 == 'electric' else '*'
        label = f"{pf1}dC({l1}),p={p1} + {pf2}dC({l2}),p={p2}"

        log_msg(log, f"Trying pair: {label}")
        try:
            results = experiment_pair(legs1, p1, t1, legs2, p2, t2)
        except Exception as e:
            results = [{'status': 'error', 'error': str(e)}]

        log['tried'].append({'exp': exp, 'result': results})

        # Check for candidates
        for r in results:
            if isinstance(r, dict) and r.get('status') == 'candidate':
                log_msg(log, f"  [CANDIDATE] c1^2={r['c1sq']}, c2^2={r['c2sq']}")

                # Full verification
                c1 = sp.sqrt(sp.S(r['c1sq']))
                c2 = sp.sqrt(sp.S(r['c2sq']))
                log_msg(log, f"  Verifying with c1={c1}, c2={c2} ...")
                ok = verify_pair(legs1, p1, t1, c1, legs2, p2, t2, c2)
                if ok:
                    sol = {
                        'form1': f"{pf1}dC({l1}), p={p1}, c={str(c1)}",
                        'form2': f"{pf2}dC({l2}), p={p2}, c={str(c2)}",
                    }
                    log['solution'] = sol
                    log_msg(log, f"  >>> SOLUTION FOUND: {sol}")
                    save_log(log)
                    return
                else:
                    log_msg(log, f"  Full verification failed (cross terms?)")
            else:
                status = r.get('status', '?') if isinstance(r, dict) else '?'
                log_msg(log, f"  -> {status}")

        save_log(log)

    log_msg(log, "--- Search complete, no solution found ---")
    save_log(log)

if __name__ == '__main__':
    run_agent()
