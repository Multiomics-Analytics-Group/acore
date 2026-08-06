# using pydantic since these are just dictionary outputs no df
from collections.abc import Callable
from typing import Any

from pydantic import BaseModel, Field


class PermutationResult(BaseModel):
    """Schema for the output of a permutation test."""

    metric: str | Callable | None = Field(
        default=None, description="Name of the metric used in the permutation test"
    )
    observed_statistic: Any = Field(
        default=None,
        description="Observed value of the metric, test statistic or scipy result object",
    )
    p_value: float = Field(description="p-value from the permutation test")
