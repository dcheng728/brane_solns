"""Index types and abstract tensor indices."""


class IndexType:
    """A vector space or tangent bundle that indices belong to.

    Parameters
    ----------
    name : str
        Human-readable name, e.g. 'Lorentz', 'Torus'.
    dim : int or sympy.Expr
        Dimension (can be symbolic).
    dummy_prefix : str
        Prefix for auto-generated dummy names (default: first letter of name).
    metric_symmetry : int
        1 = symmetric metric (default), -1 = antisymmetric, 0 = no metric.
    """

    def __init__(self, name, dim=None, dummy_prefix=None, metric_symmetry=1):
        self.name = name
        self.dim = dim
        self.dummy_prefix = dummy_prefix or name[0]
        self.metric_symmetry = metric_symmetry

    def __repr__(self):
        return f"IndexType({self.name!r}, dim={self.dim})"

    def __eq__(self, other):
        return self is other

    def __hash__(self):
        return id(self)


class Index:
    """An abstract tensor index.

    Parameters
    ----------
    name : str
        Index name, e.g. 'a', 'mu'.
    index_type : IndexType
        The vector space this index lives in.
    is_up : bool
        True = contravariant (upper), False = covariant (lower).

    Usage
    -----
    >>> L = IndexType('Lorentz', 4)
    >>> a = Index('a', L)           # upper index
    >>> -a                           # lower index (same name, flipped)
    >>> R(a, b, -c, -d)             # mixed tensor
    """

    def __init__(self, name, index_type, is_up=True):
        self.name = name
        self.index_type = index_type
        self.is_up = is_up

    def __neg__(self):
        """Flip variance: upper <-> lower."""
        return Index(self.name, self.index_type, not self.is_up)

    def __repr__(self):
        sign = '' if self.is_up else '-'
        return f"{sign}{self.name}"

    def __eq__(self, other):
        if not isinstance(other, Index):
            return NotImplemented
        return (self.name == other.name
                and self.index_type is other.index_type
                and self.is_up == other.is_up)

    def __hash__(self):
        return hash((self.name, id(self.index_type), self.is_up))

    @property
    def is_contravariant(self):
        return self.is_up

    @property
    def is_covariant(self):
        return not self.is_up

    def matches(self, other):
        """Check if this index can contract with other (same name+type, opposite variance)."""
        return (self.name == other.name
                and self.index_type is other.index_type
                and self.is_up != other.is_up)


def indices(names, index_type, is_up=True):
    """Create multiple indices at once.

    Parameters
    ----------
    names : str
        Space-separated index names, e.g. 'a b c d'.
    index_type : IndexType
    is_up : bool
        Default variance for all created indices.

    Returns
    -------
    tuple of Index

    Usage
    -----
    >>> a, b, c, d = indices('a b c d', Lorentz)
    """
    name_list = names.split()
    result = tuple(Index(n, index_type, is_up) for n in name_list)
    if len(result) == 1:
        return result[0]
    return result
