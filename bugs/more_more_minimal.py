from sympy import Array, diag, symbols
from sympy.tensor.tensor import TensorIndexType, TensorHead, tensor_indices, tensor_heads
import sympy as sp

L = TensorIndexType("L")
L2 = TensorIndexType("L2", dim=2)
i, j, k, l, m = tensor_indices("i j k l m", L)
A, B, C, D = tensor_heads("A B C D", [L])
H = TensorHead("H", [L, L])
M = TensorHead("M", [L]*3)
K = TensorHead("K", [L]*4)

# # Tensors with contractions in replacements:
# expr = K(i, j, k, -k)
# repl = {K(i, j, k, -k): [[1, 1], [1, 1]]}
# assert expr._extract_data(repl) == ([i, j], Array([[1, 1], [1, 1]]))

# expr = K(i, j, k, -k)
# repl = {K(i, j, l,-k): sp.MutableDenseNDimArray.zeros(2, 2, 2, 2), L:diag(1,1)}
# assert expr._extract_data(repl) == ([i, j], sp.MutableDenseNDimArray.zeros(2, 2))




# Riemann tensor expression in terms of Gamma, only kept one term for minimal example
Riemann_expr = (
  + M(l, -m, -k) * M(k, -i, -j)
)

# Use Christoffel symbols for flat spacetime, which is just zeros
christoffel_array = sp.Array.zeros(2,2,2)

repl = {
    M(i, -j, -k): christoffel_array, L: sp.diag(1,1)
}

print(Riemann_expr.replace_with_arrays(repl)[0,0,0,0])

# repl = {
#     M(k, -i, -j): christoffel_array, L: sp.diag(1,1)
# }

# print(Riemann_expr.replace_with_arrays(repl)[0,0,0,0])

# repl = {
#     M(j, -k, -i): christoffel_array, L: sp.diag(1,1)
# }
# print(Riemann_expr.replace_with_arrays(repl)[0,0,0,0])