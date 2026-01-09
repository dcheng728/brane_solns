import sympy as sp

def impose_harmonic_condition_sym(expr, dimension, Hpp, Hp, radius_symbol):
    """
    Impose harmonic condition at a symbol level (algebraic)

    Should be applied after derivatives have been computed and substituted for H', H''
    """
    harmonic_sub = {
        Hpp: -(dimension-1)/radius_symbol*Hp,
    }
    return expr.subs(harmonic_sub)

def minkowski_metric(dimensions):
    return sp.diag(*([-1] + [1]*(dimensions-1)))

def euclidean_metric(dimensions):
    return sp.diag(*([1]*dimensions))

def ricci_tensor(metric, coordinates):
    n_dims = metric.shape[0]
    inv_metric = metric.inv()

    # compute christoffel symbols Order: Gamma[M,N,P] = Gamma^M{}_{NP}
    Gamma = [[[sp.S(0) for _ in range(n_dims)] for _ in range(n_dims)] for _ in range(n_dims)]
    for M in range(n_dims):
        for N in range(n_dims):
            for P in range(n_dims):
                sum_term = sp.S(0)
                for Q in range(n_dims):
                    d_gNQ_dXP = sp.diff(metric[N,Q], coordinates[P])
                    d_gPQ_dXN = sp.diff(metric[P,Q], coordinates[N])
                    d_gNP_dXQ = sp.diff(metric[N,P], coordinates[Q])
                    sum_term += inv_metric[M,Q] * (d_gNQ_dXP + d_gPQ_dXN - d_gNP_dXQ)
                Gamma[M][N][P] = sp.simplify(sp.Rational(1,2) * sum_term)
        print(f"Computed Gamma^{M}_...", end=",")

    # compute ricci tensor
    # R_{MN} = \partial_P \Gamma^P_{MN} - \partial_{N}\Gamma^P_{MP} + \Gamma^P_{MN}\Gamma^Q_{PQ} - \Gamma^P_{MQ}\Gamma^Q_{NP}
    Ricci_MN = sp.zeros(n_dims,n_dims)
    for M in range(n_dims):
        for N in range(n_dims):
            sum_term = sp.S(0)
            # derivative terms
            for P in range(n_dims):
                d_GammaPMN_dXP = sp.diff(Gamma[P][M][N], coordinates[P])
                d_GammaPMP_dXN = sp.diff(Gamma[P][M][P], coordinates[N])
                sum_term += d_GammaPMN_dXP - d_GammaPMP_dXN

            # quadratic terms
            for P in range(n_dims):
                for Q in range(n_dims):
                    sum_term += Gamma[P][M][N] * Gamma[Q][P][Q] - Gamma[P][M][Q] * Gamma[Q][N][P]

            Ricci_MN[M,N] = sum_term
        print(f"Computed R[{M},...]", end=",")

    return Ricci_MN

def ricci_scalar(ricci_tensor, inv_metric):
    n_dims = ricci_tensor.shape[0]
    R = sp.S(0)
    for M in range(n_dims):
        for N in range(n_dims):
            R += inv_metric[M,N] * ricci_tensor[M,N]
    return sp.simplify(R)