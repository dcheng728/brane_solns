# working with type IIB effective actions in sympy

Because sympy generically does not have support for working with tensors (it does but that module does not suit our task here very well).
The philosophy is to work with plain sympy function and variables.
Because covariant derivatives are distributive and satisfy the Leibniz rule, when it comes to working with their algebra we can treat them as plain derivatives.
So the approach is to 
1. derive expressions by hand
2. put all derivatives in covariant form
3. translate them into sympy as plain derivatives
4. perform all covariant-derivative-related operations we need
5. translate the plain derivatives to covariant derivatives.

For example, in sympy, we can define the indices as derivatives

```python 
import sympy as sp

m,n,r,s = sp.symbols('m n r s', real=True)
tau1 = sp.Function('tau_1', real=True)(m, n, r, s)
tau2 = sp.Function('tau_2', real=True, positive=True)(m, n, r, s)
tau = tau1 + sp.I*tau2
taub = tau1 - sp.I*tau2
```

then all the distributive and leibniz properties work the same way!