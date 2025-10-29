import numpy as np
import itertools
from collections.abc import Iterable


def jaccard_similarity(set1: set, set2: set) -> float:
    """
    Compute the Jaccard similarity between two sets.

    Parameters
    ----------
    set1 : set
        First set of nodes.
    set2 : set
        Second set of nodes.
    Returns
    -------
    float
        Jaccard similarity coefficient,
        which is the size of the intersection divided by
        the size of the union of the two sets.

    Example
    -------
    >>> set1 = {1, 2, 3}
    >>> set2 = {2, 3, 4}
    >>> jaccard_similarity(set1, set2)
    0.5
    """
    intersection = len(set1 & set2)
    union = len(set1 | set2)

    if union != 0:
        return intersection / union
    else:
        return 0


def avg_jaccard(group: Iterable) -> tuple:
    """
     Compute the average pairwise Jaccard similarity within a group of graphs.

     Parameters
     ----------
     group : Iterable
         An iterator of iterable objects, such as a list of sets.

    Returns
     -------
     tuple
         A tuple containing the average Jaccard similarity and its standard deviation.

     Example
     -------
     >>> group = [{1, 2, 3}, {2, 3, 4}, {3, 4, 5}]
     >>> avg_jaccard(group)
     (0.3333333333333333, 0.0)
    """

    ## PRECONDITIONS
    if not isinstance(group, Iterable):
        raise TypeError("Input must be an iterable of iterables.")
    if len(group) > 0:
        for g in group:
            if not isinstance(g, Iterable):
                raise TypeError(
                    "Each element in the group must be an iterable (set, list, or tuple)."
                )
    elif len(group) == 0:
        raise ValueError("Input group is empty. Cannot compute Jaccard similarity.")

    # init the list
    scores = []

    for g1, g2 in itertools.combinations(group, 2):
        # get sets
        set1 = set(g1)
        set2 = set(g2)
        # compute Jaccard similarity
        scores.append(jaccard_similarity(set1, set2))

    if scores:
        # if scores is not empty, return the average
        return np.mean(scores), np.std(scores)
    # if scores is empty (less than 2 in group), return 0
    else:
        return 0, 0


def btwn_jaccard(group1: Iterable, group2: Iterable) -> tuple:
    """
    Compute the average Jaccard similarity between two groups of graphs.

    Parameters
    ----------
    group1 : Iterable
        An iterator of iterable objects, such as a list of sets.
    group2 : Iterable
        Another iterator of iterable objects.

    Returns
    -------
    tuple
        A tuple containing the average Jaccard similarity and its standard deviation.

    Example
    -------
    >>> group1 = [{1, 2, 3}, {2, 3, 4}]
    >>> group2 = [{3, 4, 5}, {4, 5, 6}]
    >>> btwn_jaccard(group1, group2)
    (0.3333333333333333, 0.0)

    """

    if not isinstance(group1, Iterable) or not isinstance(group2, Iterable):
        raise TypeError("Both inputs must be iterables of iterables.")
    if len(group1) == 0 or len(group2) == 0:
        raise ValueError("Both input groups must be non-empty.")

    # init the list
    scores = []

    for g1, g2 in itertools.product(group1, group2):
        # check
        if not isinstance(g1, Iterable) or not isinstance(g2, Iterable):
            raise TypeError(
                "Each element in the groups must be an iterable (set, list, or tuple)."
            )

        # get sets
        set1 = set(g1)
        set2 = set(g2)

        # compute Jaccard similarity
        scores.append(jaccard_similarity(set1, set2))

    if scores:
        # if scores is not empty, return the average
        return np.mean(scores), np.std(scores)
    # if scores is empty (less than 2 in group), return 0
    else:
        return 0, 0
