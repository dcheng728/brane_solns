from typing import List
from tensorGlow.core.index import Index, IndexType

class Expr:
    """Base class for all tensor expressions."""

class IndexedObject:
    def __init__(self, 
                 name: str,
                 indices: List[Index]):
        self.name = name
        self.index_types = [index.index_type for index in indices]
        self.indices = indices
        self.rank = len(indices)

class Tensor(IndexedObject):
    