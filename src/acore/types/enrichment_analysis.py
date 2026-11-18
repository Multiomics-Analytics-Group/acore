from pandera.pandas import DataFrameModel, Field


class EnrichmentAnalysisSchema(DataFrameModel):
    """
    Schema for the enrichment analysis results DataFrame.
    """

    terms: str = Field(nullable=False)
    identifiers: str = Field(nullable=False)
    foreground: int = Field(ge=0)
    background: int = Field(ge=0)
    foreground_pop: int = Field(gt=0)
    background_pop: int = Field(gt=0)
    odds_ratio: float = Field(ge=0, nullable=True)
    odds_ratio_lower: float = Field(ge=0, nullable=True)
    pvalue: float = Field(ge=0, le=1)
    padj: float = Field(ge=0, le=1)
    rejected: bool = Field(nullable=False)
    direction: str
    comparison: str
