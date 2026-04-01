import sympy as sp

class IndexType:
    def __init__(self, name, dimension, dummy_index ):
        self.name = name
        self.dimension = dimension
        self.dummy_index = dummy_index


    def __repr__(self):
        return f"IndexType(name={self.name}, dimension={self.dimension}, dummy_index={self.dummy_index})"
    
class Index:
    def __init__(self, 
                 name, 
                 index_type: IndexType,
                 up: bool = True):
        self.name = name
        self.index_type = index_type
        self.up = up

    def __neg__(self):
        return Index(self.name, self.index_type, not self.up)

    def __repr__(self):
        return f"Index(name={self.name}, index_type={self.index_type})"