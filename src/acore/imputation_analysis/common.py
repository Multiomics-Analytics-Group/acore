import logging
from collections.abc import Iterable

import pandas as pd

logger = logging.getLogger(__name__)


def select_imputation_targets(
    data: pd.DataFrame,
    on_cols: Iterable[str] | None,
    on_rows: Iterable[str] | None,
    drop_cols: Iterable[str] | None,
) -> tuple[pd.DataFrame, list[str], list[str] | None]:
    """Copy `data`, drop `drop_cols` and resolve which columns and rows to impute.

    :return: tuple of the copied DataFrame (without `drop_cols`), the numeric columns
                     to impute and the row labels found in the index (`None` if `on_rows`
                     is `None`).
    :raises TypeError: if `on_cols` contains non-numeric columns.
    """
    if drop_cols:
        if on_cols is not None:
            overlap = set(on_cols) & set(drop_cols)
            if overlap:
                logger.warning(
                    f"Columns in both on_cols and drop_cols will be dropped, not filled: {overlap}"
                )
        df = data.drop(columns=drop_cols)  # already a copy
    else:
        df = data.copy()

    if on_cols is None:
        for col in df.columns[df.dtypes == object]:
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError):
                pass
        cols = df.select_dtypes(include="number").columns.tolist()
        non_numeric = df.select_dtypes(exclude="number").columns.tolist()
        if non_numeric:
            logger.warning(f"Non-numeric columns ignored for imputation: {non_numeric}")
    else:
        cols = [c for c in on_cols if c in df.columns]
        non_numeric = [c for c in cols if not pd.api.types.is_numeric_dtype(df[c])]
        if non_numeric:
            raise TypeError(f"Non-numeric columns passed to `on_cols`: {non_numeric}")

    rows = None
    if on_rows is not None:
        on_rows = list(on_rows)
        rows = [r for r in on_rows if r in df.index]
        if len(rows) != len(on_rows):
            logger.warning(
                f"Some rows in `on_rows` were not found in the DataFrame index and will be skipped: "
                f"{set(on_rows) - set(rows)}"
            )
    return df, cols, rows
