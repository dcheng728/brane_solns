"""Standardized display for problem sheets and derivations."""


class Problem:
    """A problem with optional derivation steps and a final answer.

    Usage
    -----
    >>> with Problem("Q3(a)", "Christoffel symbol") as p:
    ...     p.step("Defining formula for the connection",
    ...            "Gamma^c_{ab}", geom.christoffel_formula(c, -a, -b))
    ...     p.step("Symmetry in lower indices",
    ...            "Gamma^c_{ba}", geom.christoffel_formula(c, -b, -a))
    ...     p.answer("Gamma^c_{ab}", geom.christoffel_formula(c, -a, -b))
    """

    _width = 60

    def __init__(self, label, title=""):
        self.label = label
        self.title = title
        self._steps = []
        self._answer_lhs = None
        self._answer_rhs = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._render()
        return False

    def step(self, description, lhs=None, rhs=None):
        """Record a derivation step.

        Parameters
        ----------
        description : str
            What this step does.
        lhs : str, optional
            Left-hand side label (e.g. "Gamma^c_{ab}").
        rhs : TensorExpr or any, optional
            The computed expression.
        """
        self._steps.append((description, lhs, rhs))

    def answer(self, lhs, rhs):
        """Record the final answer.

        Parameters
        ----------
        lhs : str
            Left-hand side label.
        rhs : TensorExpr or any
            The computed result.
        """
        self._answer_lhs = lhs
        self._answer_rhs = rhs

    def _render(self):
        # Header
        header = f"{self.label}: {self.title}" if self.title else self.label
        print()
        print("=" * self._width)
        print(header)
        print("=" * self._width)

        # Steps
        for i, (desc, lhs, rhs) in enumerate(self._steps, 1):
            print(f"\n  Step {i}: {desc}")
            if lhs is not None and rhs is not None:
                print(f"    {lhs} = {rhs}")
            elif rhs is not None:
                print(f"    {rhs}")
            elif lhs is not None:
                print(f"    {lhs}")

        # Answer
        if self._answer_lhs is not None:
            if self._steps:
                print()
                print(f"  Result:")
            else:
                print()
            print(f"    {self._answer_lhs} = {self._answer_rhs}")
