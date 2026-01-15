from sympy.tensor.tensor import TensorIndexType, TensorHead, tensor_indices # Tested with sympy 1.14.0
import sympy as sp

# Define tensor indices
Lorentz = TensorIndexType("Lorentz", dummy_name=r"\lambda")
Lorentz.set_metric(sp.diag(-1,1,1,1))
mu, nu, rho, alpha, beta = tensor_indices(r"\mu \nu \rho \alpha \beta", Lorentz)

# Define Christoffel symbols and Riemann tensor
Gamma = TensorHead(r"\Gamma", [Lorentz,Lorentz,Lorentz]) # Christoffel symbols
Riemann_tensor = TensorHead(r"R", [Lorentz,Lorentz,Lorentz,Lorentz]) # Riemann tensor 

# Riemann tensor expression in terms of Gamma, only kept one term for minimal example
Riemann_expr = (
  + Gamma(alpha, -beta, -rho) * Gamma(rho, -mu, -nu)
)

# Use Christoffel symbols for flat spacetime, which is just zeros
christoffel_array = sp.MutableDenseNDimArray.zeros(4,4,4)

# If you replace tensor like this, it works
repl = {
    Gamma(mu,-nu,-rho): christoffel_array
}
print(Riemann_expr.replace_with_arrays(repl)[0,0,0,0])

# But if you replace with indices in a different order
# like this, it gives ValueError: incompatible indices: [\lambda_0, -\mu, -\nu] and [\lambda_0, -\rho, -\mu]
repl = {
    Gamma(nu,-rho, -mu): christoffel_array
}

print(Riemann_expr.replace_with_arrays(repl)[0,0,0,0])
