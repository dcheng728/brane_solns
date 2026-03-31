"""Extensions: dimensional reduction, block metrics."""

from .dimensional_reduction import IndexSplit
from .block_metric import BlockMetric, ScalarField

__all__ = [
    'IndexSplit',
    'BlockMetric', 'ScalarField',
]
