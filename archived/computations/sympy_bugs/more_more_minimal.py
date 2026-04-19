from sympy.tensor.array import Array
from sympy.tensor.tensor import TensorIndexType, TensorHead, tensor_indices, tensor_heads
from sympy.matrices import diag

L = TensorIndexType("L")
i, j, k, l, m = tensor_indices("i j k l m", L)
M = TensorHead("M", [L]*3)

# Product of rank-3 tensors with contractions in replacements
expr = M(l, -m, -k) * M(k, -i, -j) 
# repl = {M(i, -k, -l): Array.zeros(2,2,2), L: diag(1, 1)} # doesn't work
repl = {M(l, -i, -j): Array.zeros(2,2,2), L: diag(1, 1)}
assert expr._extract_data(repl) == ([l, -m, -i, -j], Array.zeros(2,2,2,2))