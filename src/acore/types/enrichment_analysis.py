from pandera.pandas import DataFrameModel, Field


class EnrichmentAnalysisSchema(DataFrameModel):

    terms: str = Field(nullable=False)
    identifiers: str = Field(nullable=False)
    foreground: int = Field(gt=0)
    background: int = Field(gt=0)
    foreground_pop: int = Field(gt=0)
    background_pop: int = Field(gt=0)
    pvalue: float = Field(gt=0, le=1)
    padj: float = Field(gt=0, le=1)
    rejected: bool = Field(nullable=False)
    direction: str
    comparison: str
