"""IO module: console and LaTeX rendering for tensor expressions.

Handles all expression types: Scalar, Tensor, Apply, Prod, Sum.
"""

import sympy as sp
from ..core.expr import Expr, Scalar, Tensor, Apply, Prod, Sum


# ═════════════════════════════════════════════════════════════════════
# Console rendering (plain text with ^{} and _{})
# ═════════════════════════════════════════════════════════════════════

def to_console(expr):
    """Render an expression for console output.

    Uses ^{ab} and _{cd} notation for upper/lower indices.
    """
    return _render(expr, _ConsoleFormatter())


# ═════════════════════════════════════════════════════════════════════
# LaTeX rendering
# ═════════════════════════════════════════════════════════════════════

def to_latex(expr, index_names=None):
    """Render an expression as LaTeX.

    Parameters
    ----------
    index_names : dict, optional
        Map index name -> LaTeX name, e.g. {'mu': r'\\mu'}.
    """
    return _render(expr, _LaTeXFormatter(index_names or {}))


# ═════════════════════════════════════════════════════════════════════
# Renderer dispatch
# ═════════════════════════════════════════════════════════════════════

def _render(expr, fmt):
    """Dispatch rendering by expression type."""
    if isinstance(expr, Sum):
        return _render_sum(expr, fmt)
    if isinstance(expr, Prod):
        return _render_prod(expr, fmt)
    if isinstance(expr, Apply):
        return _render_apply(expr, fmt)
    if isinstance(expr, Tensor):
        return _render_tensor(expr, fmt)
    if isinstance(expr, Scalar):
        return fmt.scalar(expr.expr)
    return repr(expr)


def _render_tensor(tensor, fmt):
    """Render a named tensor with grouped upper/lower indices."""
    return fmt.indexed(tensor.head.name, tensor.indices)


def _render_apply(apply_expr, fmt):
    """Render operator application: partial_{a}(g_{b c})."""
    op_name = apply_expr.op.name
    idx = apply_expr.deriv_index
    operand_str = _render(apply_expr.operand, fmt)
    return fmt.apply(op_name, idx, operand_str)


def _render_prod(prod, fmt):
    """Render a product of factors with a coefficient."""
    if not prod.factors:
        return fmt.scalar(prod.coeff)

    parts = [_render(f, fmt) for f in prod.factors]
    factors_str = fmt.product_join(parts)

    if prod.coeff == 1:
        return factors_str
    if prod.coeff == -1:
        return fmt.negate(factors_str)
    return fmt.scaled(prod.coeff, factors_str)


def _render_sum(sum_expr, fmt):
    """Render a sum of terms."""
    if not sum_expr.terms:
        return "0"

    parts = []
    for i, term in enumerate(sum_expr.terms):
        term_str = _render(term, fmt)
        if i == 0:
            parts.append(term_str)
        elif term_str.startswith('-'):
            parts.append(f" {term_str}")
        else:
            parts.append(f" + {term_str}")

    return ''.join(parts)


# ═════════════════════════════════════════════════════════════════════
# Formatter classes
# ═════════════════════════════════════════════════════════════════════

class _ConsoleFormatter:
    """Plain console formatter: R(a, b, -c, -d), partial(-a)(g(-b, -c))."""

    def indexed(self, name, indices):
        if not indices:
            return name
        idx_strs = []
        for idx in indices:
            if idx.is_up:
                idx_strs.append(idx.name)
            else:
                idx_strs.append(f'-{idx.name}')
        return f"{name}({', '.join(idx_strs)})"

    def apply(self, op_name, idx, operand_str):
        if idx.is_up:
            return f"{op_name}({idx.name})({operand_str})"
        else:
            return f"{op_name}(-{idx.name})({operand_str})"

    def product_join(self, parts):
        return ' '.join(parts)

    def scalar(self, expr):
        return str(expr)

    def negate(self, s):
        return f"-{s}"

    def scaled(self, coeff, s):
        return f"{coeff} {s}"


class _LaTeXFormatter:
    """LaTeX formatter: R^{a b}_{c d}."""

    def __init__(self, index_names=None):
        self.index_names = index_names or {}

    def _idx_name(self, idx):
        return self.index_names.get(idx.name, idx.name)

    def indexed(self, name, indices):
        if not indices:
            return name
        groups = _group_indices(indices, self.index_names)
        # Escape underscores in tensor names for LaTeX
        latex_name = name.replace('_', r'\_') if '_' in name and '\\' not in name else name
        result = latex_name
        for is_up, names in groups:
            joined = ' '.join(names)
            if is_up:
                result += f'^{{{joined}}}'
            else:
                result += f'_{{{joined}}}'
        return result

    def apply(self, op_name, idx, operand_str):
        idx_name = self._idx_name(idx)
        if op_name == 'partial':
            return rf'\partial_{{{idx_name}}}\left({operand_str}\right)'
        return rf'{op_name}_{{{idx_name}}}\left({operand_str}\right)'

    def product_join(self, parts):
        return r'\,'.join(parts)

    def scalar(self, expr):
        return sp.latex(expr)

    def negate(self, s):
        return f"-{s}"

    def scaled(self, coeff, s):
        return rf'{sp.latex(coeff)}\,{s}'


# ═════════════════════════════════════════════════════════════════════
# Helpers
# ═════════════════════════════════════════════════════════════════════

def _group_indices(indices, name_map=None):
    """Group consecutive upper/lower indices.

    Returns list of (is_up, [name_strings]).
    """
    if name_map is None:
        name_map = {}
    groups = []
    for idx in indices:
        name = name_map.get(idx.name, idx.name)
        if groups and groups[-1][0] == idx.is_up:
            groups[-1][1].append(name)
        else:
            groups.append((idx.is_up, [name]))
    return groups
