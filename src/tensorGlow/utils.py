"""Shared utilities for tensorGlow."""

import itertools


_dummy_counters = {}


def fresh_dummy_name(prefix='L'):
    """Generate a fresh dummy index name like L_0, L_1, ..."""
    count = _dummy_counters.get(prefix, 0)
    _dummy_counters[prefix] = count + 1
    return f"{prefix}_{count}"


def reset_dummy_counters():
    """Reset all dummy name counters (useful in tests)."""
    _dummy_counters.clear()


def perm_sign(perm):
    """Sign of a permutation given as a list.

    Parameters
    ----------
    perm : list of int
        A permutation of [0, 1, ..., n-1].

    Returns
    -------
    int : +1 or -1
    """
    n = len(perm)
    visited = [False] * n
    sign = 1
    for i in range(n):
        if visited[i]:
            continue
        j = i
        cycle_len = 0
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        if cycle_len % 2 == 0:
            sign *= -1
    return sign
