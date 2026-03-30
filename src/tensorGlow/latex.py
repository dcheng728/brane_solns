"""LaTeX output for tensor expressions."""

from .expr import TensorAtom, TensorProduct, TensorSum, ScalarExpr
import sympy as sp


def to_latex(expr, index_names=None):
    """Convert a TensorExpr to a LaTeX string.

    Parameters
    ----------
    expr : TensorExpr
    index_names : dict, optional
        Custom LaTeX names for indices, e.g. {'mu': r'\\mu'}.

    Returns
    -------
    str
    """
    if index_names is None:
        index_names = {}

    if isinstance(expr, TensorSum):
        return _latex_sum(expr, index_names)
    if isinstance(expr, TensorProduct):
        return _latex_product(expr, index_names)
    if isinstance(expr, ScalarExpr):
        return sp.latex(expr.expr)
    raise TypeError(f"Cannot convert {type(expr)} to LaTeX")


def _latex_index(idx, names):
    """Render a single index."""
    name = names.get(idx.name, idx.name)
    return name


def _latex_atom(atom, names):
    """Render a TensorAtom like R^{ab}_{cd}."""
    head_name = atom.head.name

    # Group consecutive upper/lower indices
    groups = []  # list of (is_up, [index_strs])
    for idx in atom.indices:
        latex_name = _latex_index(idx, names)
        if groups and groups[-1][0] == idx.is_up:
            groups[-1][1].append(latex_name)
        else:
            groups.append((idx.is_up, [latex_name]))

    result = head_name
    for is_up, idx_strs in groups:
        joined = ' '.join(idx_strs)
        if is_up:
            result += f"^{{{joined}}}"
        else:
            result += f"_{{{joined}}}"

    return result


def _latex_product(product, names):
    """Render a TensorProduct."""
    if not product.atoms:
        return sp.latex(product.coeff)

    atoms_str = r'\,'.join(_latex_atom(a, names) for a in product.atoms)

    if product.coeff == 1:
        return atoms_str
    if product.coeff == -1:
        return f"-{atoms_str}"

    coeff_str = sp.latex(product.coeff)
    return f"{coeff_str}\\,{atoms_str}"


def _latex_sum(tensor_sum, names):
    """Render a TensorSum."""
    if not tensor_sum.terms:
        return "0"

    parts = []
    for i, term in enumerate(tensor_sum.terms):
        term_str = _latex_product(term, names)
        if i > 0 and not term_str.startswith('-'):
            parts.append(f" + {term_str}")
        elif i > 0:
            parts.append(f" {term_str}")
        else:
            parts.append(term_str)

    return ''.join(parts)
