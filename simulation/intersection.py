"""
Intersection computation for pairwise carrier overlap.

O(n + m) two-pointer algorithm.
"""

from __future__ import annotations


def intersection_size(a: list[int], b: list[int]) -> int:
    """Two-pointer O(n + m) intersection size. Arrays must be sorted."""
    i, j = 0, 0
    count = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            count += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    return count


def is_eligible(a: list[int], b: list[int]) -> bool:
    """Eligible = no shared pathogenic loci (empty intersection)."""
    return intersection_size(a, b) == 0
