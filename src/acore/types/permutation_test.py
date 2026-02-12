# using pydantic since these are just dictionary outputs no df
from typing import Any, Optional, Callable, Union
from pydantic import BaseModel, Field


class PermutationResult(BaseModel):
    metric: Optional[Union[str, Callable]] = Field(
        default=None, description="Name of the metric used in the permutation test"
    )
    observed: Any = Field(default=None, description="Observed value of the metric")
    p_value: float = Field(description="p-value from the permutation test")
