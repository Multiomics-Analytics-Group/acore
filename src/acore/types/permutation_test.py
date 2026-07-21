# using pydantic since these are just dictionary outputs no df
from typing import Any, Callable, Optional, Union

from pydantic import BaseModel, Field


class PermutationResult(BaseModel):
    """Schema for the output of a permutation test."""

    metric: Optional[Union[str, Callable]] = Field(
        default=None, description="Name of the metric used in the permutation test"
    )
    observed: Any = Field(
        default=None,
        description="Observed value of the metric, test statistic or scipy result object",
    )
    p_value: float = Field(description="p-value from the permutation test")
