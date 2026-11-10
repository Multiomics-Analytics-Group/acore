import pandera.pandas as pa
import pydantic
from pandera.pandas import DataFrameModel, Field, SeriesSchema

# from pandera.typing.pandas import DataFrame, Series


class TwoComponentSchema(DataFrameModel):
    """
    Schema for the PCA components DataFrame.
    """

    group: str | None = Field(nullable=False)
    x: float = Field(coerce=True)
    y: float = Field(coerce=True)


class TwoLoadingsSchema(TwoComponentSchema):
    """
    Schema for the PCA loadings DataFrame.
    """

    value: float | None = Field(
        nullable=False,
        description="Variance of feature explained by the extracted components.",
    )


class ThreeComponentSchema(TwoComponentSchema):
    """
    Schema for the components DataFrame if three (or more) components are requested.
    Further components are named "3", "4", ... and are not validated.
    """

    z: float = Field(coerce=True)

    class Config:
        strict = False  # allow additional components and annotation columns


class ThreeLoadingsSchema(ThreeComponentSchema):
    """
    Schema for the PCA loadings DataFrame if three (or more) components are requested.
    """

    value: float | None = Field(
        nullable=False,
        description="Variance of feature explained by the extracted components.",
    )


TwoVariance = SeriesSchema(float, checks=pa.Check(lambda x: x >= 0, element_wise=True))

# explained variance ratios are fractions of the total variance
VarianceRatio = SeriesSchema(
    float, checks=pa.Check(lambda x: 0 <= x <= 1, element_wise=True)
)

AnnotationSchema = SeriesSchema(str)


class AnnotationResult(pydantic.BaseModel):
    """
    Represents the annotation result from exploratory analysis.
    """

    x_title: str
    y_title: str
    z_title: str | None = pydantic.Field(default=None)
    group: str | None = pydantic.Field(default=None)
