from pandera.pandas import DataFrameModel, Field


class AnovaSchema(DataFrameModel):
    """
    Schema for the enrichment analysis results DataFrame.
    """

    group1: str = Field(nullable=False)
    group2: str = Field(nullable=False)
    mean_group1: float = Field(alias="mean(group1)", nullable=True)
    std_group1: float = Field(ge=0, alias="std(group1)", nullable=True)
    mean_group2: float = Field(alias="mean(group2)", nullable=True)
    std_group2: float = Field(ge=0, alias="std(group2)", nullable=True)
    # p-values
    # pvalue: float = Field(alias="posthoc pvalue")
    t_statistics: float = Field(alias="T-Statistics", nullable=True)
    pvalue: float = Field(alias="pvalue", nullable=True)

    # fold change
    log2FC: float = Field(nullable=True)
    FC: float = Field(nullable=True)
    # test statistic

    # pvalue: float = Field(ge=0, le=1)
    # correction for multiple testing
    padj: float = Field(ge=0, le=1)
    correction: str = Field()
    rejected: bool = Field(
        nullable=False
    )  # design decision to not allow nullable values, should be set to zero in case of
    # NaN?
    neg_log10_p_value: float = Field(alias="-log10 pvalue")
    Method: str = Field()


class AnovaSchemaMultiGroup(AnovaSchema):
    """
    Schema for more than two groups
    """

    t_statistics: float = Field(alias="posthoc T-Statistics")
    posthoc_pvalue: float = Field(alias="posthoc pvalue")
    f_statistics: float = Field(alias="F-statistics")
    posthoc_padj: float = Field(alias="posthoc padj", ge=0, le=1)

    posthoc_paired: bool = Field(alias="posthoc Paired", ge=0, le=1)
    posthoc_parametric: bool = Field(alias="posthoc Parametric", ge=0, le=1)
    posthoc_dof: float = Field(alias="posthoc dof", ge=0)
    posthoc_tail: str = Field(alias="posthoc tail")
    posthoc_BF10: str = Field(alias="posthoc BF10")
    posthoc_effsize: float = Field(alias="posthoc effsize")
    efftype: str = Field(alias="efftype")


class AncovaSchema(DataFrameModel):
    """
    Schema for the enrichment analysis results DataFrame.
    """

    t_statistics: float = Field(alias="posthoc T-Statistics")
    posthoc_pvalue: float = Field(alias="posthoc pvalue")
    coef: float = Field()
    std_err: float = Field(alias="std err")
    conf_int_low: float = Field(alias="Conf. Int. Low")
    conf_int_upp: float = Field(alias="Conf. Int. Upp.")
    f_statistics: float = Field(alias="F-statistics")
    posthoc_padj: float = Field(alias="posthoc padj", ge=0, le=1)
