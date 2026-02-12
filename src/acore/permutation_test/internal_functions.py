import numpy as np
import pandas as pd


def _permute(*groups, rng=np.random.default_rng(seed=12345)):
    """
    Perform a single permutation of the groups.

    Parameters
    ----------
    *groups : iterable
        The groups to permute.
    rng : np.random.Generator, optional
        The random number generator to use for shuffling.
        Defaults to a new random number generator.

    Returns
    -------
    list
        The permuted groups.
    """
    # shuffle the combined array
    combined = np.concatenate(groups)
    rng.shuffle(combined)
    # split the shuffled array into new groups of og size
    new_groups = []
    start = 0
    for group in groups:
        end = start + len(group)
        new_groups.append(combined[start:end])
        start = end
    return new_groups


def _contingency_table(*groups, to_np: bool = True) -> np.array:
    """
    Create a contingency table from the provided groups.

    Parameters
    ----------
    *groups : iterable
        The groups to include in the contingency table.
    to_np : bool, optional
        Whether to return the table as a NumPy array (default is True).

    Returns
    -------
    np.array
        The contingency table as a NumPy array or a pandas DataFrame.
    """
    # PRECONDITIONS
    for group in groups:
        if not isinstance(group, (list, pd.Series, np.ndarray)):
            raise TypeError("Each group must be a list, pandas Series, or numpy array.")

    # MAIN FUNCTION
    # Create contingency table
    table = (
        pd.concat([pd.Series(group).value_counts() for group in groups], axis=1)
        .fillna(0)
        .T
    )

    if to_np:
        return table.to_numpy()
    else:
        return table


def _check_degeneracy(array: np.ndarray) -> bool:
    """
    Check for degeneracy in the differences of paired samples.

    Parameters
    ----------
    array : np.ndarray
        The differences array. e.g. cond1 - cond2

    Returns
    -------
    bool
        True if the conditions are degenerate, False otherwise.
    """
    return np.all(array == array[0])
