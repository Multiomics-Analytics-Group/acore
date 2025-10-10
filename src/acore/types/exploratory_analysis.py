from typing import Optional

import pandera.pandas as pa
import pydantic
from pandera.pandas import DataFrameModel, Field, SeriesSchema

# from pandera.typing.pandas import DataFrame, Series


class TwoComponentSchema(DataFrameModel):
    """
    Schema for the PCA components DataFrame.
    """

    group: Optional[str] = Field(nullable=False)
    x: float = Field(coerce=True)
    y: float = Field(coerce=True)


class TwoLoadingsSchema(TwoComponentSchema):
    """
    Schema for the PCA loadings DataFrame.
    """

    value: Optional[float] = Field(
        nullable=False,
        description="Variance of feature explained by the extracted components.",
    )


TwoVariance = SeriesSchema(float, checks=pa.Check(lambda x: x >= 0, element_wise=True))

AnnotationSchema = SeriesSchema(str)

# Not used currently - integrate pandera types with pydantic
# https://pandera.readthedocs.io/en/stable/pydantic_integration.html
# class PcaResult(pydantic.BaseModel):
#     """
#     Represents the result of a Principal Component Analysis (PCA).
#     """

#     components: DataFrame[TwoComponentSchema]
#     loadings: DataFrame[TwoLoadingsSchema]
#     variance: TwoVariance


class AnnotationResult(pydantic.BaseModel):
    """
    Represents the annotation result from exploratory analysis.
    """

    x_title: str
    y_title: str
    group: Optional[str] = pydantic.Field(default=None)
