"""tensorGlow — Abstract tensor algebra for differential geometry and physics."""

from .index import IndexType, Index, indices
from .symmetry import TensorSymmetry
from .tensor_head import TensorHead
from .expr import TensorExpr, TensorAtom, TensorProduct, TensorSum, ScalarExpr
from .canonicalize import canonicalize_expr, canonicalize_product, rename_dummies
from .metric import MetricTensor
from .derivative import CovDerivative, DerivativeAtom, PartialDerivative
from .curvature import RiemannGeometry
from .dimensional_reduction import IndexSplit
from .block_metric import BlockMetric, ScalarField
from .latex import to_latex

__all__ = [
    'IndexType', 'Index', 'indices',
    'TensorSymmetry',
    'TensorHead',
    'TensorExpr', 'TensorAtom', 'TensorProduct', 'TensorSum', 'ScalarExpr',
    'canonicalize_expr', 'canonicalize_product', 'rename_dummies',
    'MetricTensor',
    'CovDerivative', 'DerivativeAtom', 'PartialDerivative',
    'RiemannGeometry',
    'IndexSplit',
    'BlockMetric', 'ScalarField',
    'to_latex',
]
