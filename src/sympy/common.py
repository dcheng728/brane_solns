import sympy as sp

from itertools import product as iterprod

# 1. Set up the fundamental variables (indices) and functions

# Define indices — change this one line to add/remove external coordinates
INDEX_LABELS = ['m', 'n', 'r', 's', 'p'] # add more if needed

indices = list(sp.symbols(' '.join(INDEX_LABELS), real=True))
index_names = dict(zip(indices, INDEX_LABELS))
# make each index available as a module-level variable (m, n, r, s, ...)
for _idx, _label in zip(indices, INDEX_LABELS):
    globals()[_label] = _idx

# Define functions (depend on all indices)
tau1 = sp.Function('tau_1', real=True)(*indices)
tau2 = sp.Function('tau_2', real=True, positive=True)(*indices)
tau = tau1 + sp.I*tau2
taub = tau1 - sp.I*tau2

H3 = sp.Function('H3', real=True)(*indices)
F3 = sp.Function('F3', real=True)(*indices)
funcs = [tau1, tau2, H3, F3]
func_names = {tau1: r'\tau_{1}', tau2: r'\tau_{2}', H3: r'H_{3}', F3: r'F_{3}'}

# SYM_Dn:  lookup dicts mapping (func, coord...) -> display Symbol
#   SYM_D0[f]            -> Symbol for f          e.g. SYM_D0[tau1] = \tau_{1}
#   SYM_D1[(f, xi)]      -> Symbol for df/dxi     e.g. SYM_D1[(tau1, m)] = \nabla_{m} \tau_{1}
#   SYM_D2[(f, (xi,xj))] -> Symbol for d^2f/dxidxj (indices sorted alphabetically)
#
# SUBS_Dn: substitution lists for use with expr.subs(), ordered highest-derivative-first:
#   expr.subs(SUBS_D2).subs(SUBS_D1).subs(SUBS_D0)
SYM_D0 = {}
SUBS_D0 = []
SYM_D1 = {}
SUBS_D1 = []
SYM_D2 = {}
SUBS_D2 = []

# Function symbols
for f in funcs:
    sym = sp.Symbol(func_names[f], real=True)
    SYM_D0[f] = sym
    SUBS_D0.append((f, sym))

# First derivative symbols
for f, xi in iterprod(funcs, indices):
    deriv = sp.diff(f, xi)
    if deriv == 0:
        continue
    name = rf'\nabla_{{{index_names[xi]}}} {func_names[f]}'
    sym = sp.Symbol(name, real=True)
    SYM_D1[(f, xi)] = sym
    SUBS_D1.append((deriv, sym))

# Second derivative symbols
for f, (xi, xj) in iterprod(funcs, iterprod(indices, indices)):
    deriv = sp.diff(f, xi, xj)
    if deriv == 0:
        continue
    key = (f, tuple(sorted([xi, xj], key=lambda x: index_names[x])))
    if key not in SYM_D2:
        ci, cj = index_names[key[1][0]], index_names[key[1][1]]
        name = rf'\nabla_{{{ci}}} \nabla_{{{cj}}} {func_names[f]}'
        SYM_D2[key] = sp.Symbol(name, real=True)
    SUBS_D2.append((deriv, SYM_D2[key]))

print(f"Function subs: {len(SUBS_D0)}")
print(f"1st deriv subs: {len(SUBS_D1)}")
print(f"2nd deriv subs: {len(SUBS_D2)}")


# Helper functions
def func2sym(expr):
    """Replace functions and their derivatives with symbols for display"""
    expr = expr.subs(SUBS_D2).subs(SUBS_D1).subs(SUBS_D0)
    return expr

def real_parts(expr):
    """Return the real part of a complex expression"""
    expr = sp.re(expr)
    expr = expr.replace(lambda e: e.func == sp.re, lambda e: e.args[0])
    expr = expr.replace(lambda e: e.func == sp.im, lambda e: sp.S.Zero)
    return expr


def deriv_monomials(*index_pairs):
    """Return derivative monomials for one or more pairs of external indices.

    Each index pair (mu, nu) contributes:
      - 2nd derivative symbols: nabla_mu nabla_nu tau_i
      - 1st derivative pairs:  (nabla_mu tau_i)(nabla_nu tau_j)

    For multiple pairs, monomials are all products across pairs.
    E.g. for pairs [(m,n), (r,s)], a monomial might be
    (nabla_m nabla_n tau_1) * (nabla_r tau_1)(nabla_s tau_2).

    Args:
        *index_pairs: tuples of (mu, nu) external index symbols

    Returns:
        (monomials, extract_coeffs) where:
          - monomials: list of sympy expressions (the monomial basis)
          - extract_coeffs(expr): returns list of tau-dependent coefficients
    """
    # For each pair, build the atomic monomials (d2 symbols and d1*d1 products)
    atoms_per_pair = []
    for mu, nu in index_pairs:
        pair_atoms = []
        # second derivatives
        for f in funcs:
            key = (f, tuple(sorted([mu, nu], key=lambda x: index_names[x])))
            if key in SYM_D2:
                pair_atoms.append(SYM_D2[key])
        # first derivative products
        d1_mu = [SYM_D1[(f, mu)] for f in funcs if (f, mu) in SYM_D1]
        d1_nu = [SYM_D1[(f, nu)] for f in funcs if (f, nu) in SYM_D1]
        for a in d1_mu:
            for b in d1_nu:
                pair_atoms.append(a * b)
        atoms_per_pair.append(pair_atoms)

    # Full monomials = cartesian product across pairs
    monomials = [sp.Mul(*combo) for combo in iterprod(*atoms_per_pair)]

    def extract_coeffs(expr):
        expr = sp.expand(expr)
        coeffs = []
        for mono in monomials:
            c = expr
            for factor in mono.as_ordered_factors():
                c = c.coeff(factor)
            coeffs.append(c)
        return coeffs

    return monomials, extract_coeffs


def decompose(target, basis, *index_pairs):
    """Decompose target as a linear combination of basis expressions.

    Solves target = sum_i c_i * basis[i] where c_i can be complex functions
    of tau_1, tau_2. All expressions should already have func2sym applied.

    Args:
        target: sympy expression to decompose
        basis: list of sympy expressions
        *index_pairs: tuples of (mu, nu) for each Riemann factor
            e.g. (m, n) for R_{ambn}, or (m, n), (r, s) for R*R products

    Returns:
        solution vector as a sympy Matrix, or empty dict if no solution.
    """
    monomials, extract_coeffs = deriv_monomials(*index_pairs)
    M = sp.Matrix([extract_coeffs(b) for b in basis]).T
    t = sp.Matrix(extract_coeffs(target))
    c = [sp.Symbol(f'c{i}') for i in range(len(basis))]
    sol = sp.solve(M @ sp.Matrix(c) - t, c)
    if not sol:
        return []
    return [sol[ci] for ci in c]


# Also define composite functions that are not independent variables
P, Pb, Q = {}, {}, {} # P_m and \bar{P}_m and Q_m
for index in indices:
    P[index] = sp.Rational(1, 2) * sp.I / tau2 * sp.diff(tau, index)
    Pb[index] = - sp.Rational(1, 2) * sp.I / tau2 * sp.diff(taub, index)
    Q[index] = - sp.Rational(1, 2) * sp.diff(tau1, index) / tau2

DP, DPb = {}, {}  # D_m P_n for index (m, n)
for mu, nu in iterprod(indices, indices):
    DP[(mu, nu)] = sp.simplify(sp.diff(P[nu], mu) - 2 * sp.I * Q[mu] * P[nu])
    DPb[(mu, nu)] = sp.simplify(sp.diff(Pb[nu], mu) + 2 * sp.I * Q[mu] * Pb[nu])

G3  = tau2**sp.Rational(-1,2) * (F3 - (tau1 + sp.I*tau2)*H3)
G3b = tau2**sp.Rational(-1,2) * (F3 - (tau1 - sp.I*tau2)*H3)
DG3, DG3b = {}, {}  # D_m G3 and D_m G3b for index m
for index in indices:
    DG3[index] = sp.simplify(sp.diff(G3, index) - sp.I * Q[index] * G3)
    DG3b[index] = sp.simplify(sp.diff(G3b, index) + sp.I * Q[index] * G3b)

# The torus metric defined with functions
g_func = sp.Matrix([[1/tau2, tau1/tau2], [tau1/tau2, (tau1**2 + tau2**2)/tau2]])
ginv_func = sp.simplify(g_func.inv())
# torus levi-civita symbol
epsilon2_down = sp.Matrix([[0, 1], [-1, 0]])
epsilon2_up = sp.simplify(ginv_func * epsilon2_down * ginv_func)

# compute Christoffel symbols using the torus metric
Gamma_a_mb = {} # For index m this gives \Gamma^a{}_{mb} = 1/2 g^{ac} g_{bc,m}
Gamma_mab = {} # For index m this gives \Gamma_{mab} = -1/2 g_{ab,m}
for index in indices:
    Gamma_a_mb[index] = sp.simplify(sp.Rational(1,2) * ginv_func * sp.diff(g_func, index))
    Gamma_mab[index] = sp.simplify(- sp.Rational(1,2) * sp.diff(g_func, index))

F4 = sp.Matrix([H3, F3]) # The doublet 3-form components
nablaF4_down = {} # For index m, this gives (\nabla_m F4)_a
for index in indices:
    nablaF4_down[index] = sp.simplify(sp.diff(F4, index) - Gamma_a_mb[index].T * F4)


# Also define the Riemann tensors, as expressions of the symbols, because one should not need to take derivatives of the Riemann tensor
Rd_mnab = {} # For indices m,n this gives R_{mnab} as 2x2 matrix in a,b
for mu, nu in iterprod(indices, indices):
    Rd_mnab[(mu, nu)] = sp.simplify( Gamma_mab[mu] * Gamma_a_mb[nu] - Gamma_a_mb[nu].T * Gamma_mab[mu].T)

Ru_mnab = {}
for mu, nu in iterprod(indices, indices):
    Ru_mnab[(mu, nu)] = sp.simplify( ginv_func * Rd_mnab[(mu, nu)] * ginv_func )

# One can show that R_{ambn} = -g_{ac}(\nabla_n \Gamma^c{}_{mb} + \Gamma^c{}_{nd}\Gamma^d{}_{mb})
Rd_ambn = {} # Gives R_{ambn} for indices m, n
for mu, nu in iterprod(indices, indices):
    Rd_ambn[(mu, nu)] = sp.simplify( (-g_func) * (sp.diff(Gamma_a_mb[mu], nu) + Gamma_a_mb[nu] * Gamma_a_mb[mu]) )

Ru_ambn = {} # For indices m,n this gives R^a{}_m{}^b{}_n as 2x2 matrix in a,b
for mu, nu in iterprod(indices, indices):
    Ru_ambn[(mu, nu)] = sp.simplify( ginv_func * Rd_ambn[(mu, nu)] * ginv_func )

Rd_ambn_sym = {} # symbols for e.g. R_{1m1n} with key (1, m, 1, n)
Ru_ambn_sym = {} # symbols for e.g. R^2{}_m{}^1{}_n with key (2, m, 1, n)
for mu, nu in iterprod(indices, indices):
    mu_n, nu_n = index_names[mu], index_names[nu]
    for a, b in iterprod(range(2), range(2)):
        Rd_ambn_sym[(a, mu, b, nu)] = sp.Symbol(rf"R_{{{a}{mu_n}{b}{nu_n}}}")
        Ru_ambn_sym[(a, mu, b, nu)] = sp.Symbol(rf"R^{{{a}}}{{}}_{{{mu_n}}}{{}}^{{{b}}}{{}}_{{{nu_n}}}")

# Give the D_m P_n and D_m Pb_n in terms of Riemann tensor symbols
DP_with_R, DPb_with_R = {}, {} # on indices (m,n) gives D_m P_n in terms of Riemann tensor symbols
for mu, nu in iterprod(indices, indices):
    A_mu_nu = -SYM_D0[tau2] * Rd_ambn_sym[0, mu, 0, nu] + sp.Rational(1,2) * sum(func2sym(ginv_func[a,b]) * Rd_ambn_sym[(a, mu, b, nu)] for a in range(2) for b in range(2))
    B_mu_nu = SYM_D0[tau1] * Rd_ambn_sym[0, mu, 0, nu] - sp.Rational(1,2) * (Rd_ambn_sym[0, mu, 1, nu] + Rd_ambn_sym[1, mu, 0, nu])
    DP_with_R[(mu, nu)] = A_mu_nu + sp.I * B_mu_nu
    DPb_with_R[(mu, nu)] = A_mu_nu - sp.I * B_mu_nu

