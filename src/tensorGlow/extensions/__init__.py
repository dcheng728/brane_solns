"""Extensions: dimensional reduction, block metrics, LaTeX output."""

from .dimensional_reduction import IndexSplit
from .block_metric import BlockMetric, ScalarField
from .latex import to_latex

__all__ = [
    'IndexSplit',
    'BlockMetric', 'ScalarField',
    'to_latex',
]
