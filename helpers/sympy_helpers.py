import sympy as sp

def find_common_factor(expr1, expr2):
    e1 = sp.cancel(expr1)   # or eq00.lhs
    e2 = sp.cancel(expr2)   # or eq00.rhs

    n1, d1 = sp.fraction(sp.together(e1))
    n2, d2 = sp.fraction(sp.together(e2))

    gens = sorted((n1.free_symbols | n2.free_symbols | d1.free_symbols | d2.free_symbols), key=lambda s: s.name)

    g_num = sp.gcd(sp.Poly(n1, *gens), sp.Poly(n2, *gens)).as_expr()
    g_den = sp.gcd(sp.Poly(d1, *gens), sp.Poly(d2, *gens)).as_expr()

    common = sp.simplify(g_num / g_den)
    print("Common factor found:", common)
    return common
