"""tensorGlow — Abstract tensor algebra for differential geometry and physics."""

from .core import (
    IndexType, Index, indices,
    TensorSymmetry,
    TensorHead,
    TensorExpr, TensorAtom, TensorProduct, TensorSum, ScalarExpr,
    canonicalize_expr, canonicalize_product, rename_dummies,
    MetricTensor,
    CovDerivative, DerivativeAtom, PartialDerivative,
    RiemannGeometry,
)

from .extensions import (
    IndexSplit,
    BlockMetric, ScalarField,
)

from .utils import (
    to_latex,
    Problem,
)

__all__ = [
    # core
    'IndexType', 'Index', 'indices',
    'TensorSymmetry',
    'TensorHead',
    'TensorExpr', 'TensorAtom', 'TensorProduct', 'TensorSum', 'ScalarExpr',
    'canonicalize_expr', 'canonicalize_product', 'rename_dummies',
    'MetricTensor',
    'CovDerivative', 'DerivativeAtom', 'PartialDerivative',
    'RiemannGeometry',
    # extensions
    'IndexSplit',
    'BlockMetric', 'ScalarField',
    # utils
    'to_latex',
    'Problem',
]
