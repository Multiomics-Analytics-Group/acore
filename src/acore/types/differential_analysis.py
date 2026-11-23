from pandera.pandas import DataFrameModel, Field


class TtestSchema(DataFrameModel):
    """
    Schema for the output of :func:`acore.differential_regulation.run_ttest`.

    One row per feature (``identifier``) comparing ``group1`` and ``group2``.
    """

    identifier: str = Field(nullable=False)
    test: str = Field()

    group1: str = Field(nullable=False)
    group2: str = Field(nullable=False)
    mean_group1: float = Field(alias="mean(group1)", nullable=True)
    std_group1: float = Field(ge=0, alias="std(group1)", nullable=True)
    mean_group2: float = Field(alias="mean(group2)", nullable=True)
    std_group2: float = Field(ge=0, alias="std(group2)", nullable=True)
    # test statistic: T for t-test, U if non-parametric (Mann-Whitney)
    t_statistics: float = Field(alias="T-Statistics", nullable=True)
    pvalue: float = Field(ge=0, le=1, nullable=True)
    # fold change (log2FC is the mean difference if data is logged)
    log2FC: float = Field(nullable=True)
    FC: float = Field(nullable=True)
    # correction for multiple testing
    padj: float = Field(ge=0, le=1, nullable=True)
    correction: str = Field()
    rejected: bool = Field(nullable=False)
    neg_log10_p_value: float = Field(alias="-log10 pvalue", nullable=True)
    Method: str = Field()


class AnovaSchema(TtestSchema):
    """
    Schema for the two-group case of :func:`acore.differential_regulation.run_anova`.

    ``run_anova`` calls ``run_ttest`` for two groups, so the result is the same.
    The ``test`` column is only reported by the t-test, so it is optional here
    (it is absent in the multi-group results inheriting from this schema).
    """

    test: str | None = Field()


class AnovaSchemaMultiGroup(AnovaSchema):
    """
    Schema for more than two groups
    """

    posthoc_comparison: str = Field(alias="posthoc comparison", nullable=False)
    t_statistics: float = Field(alias="posthoc T-Statistics")
    posthoc_pvalue: float = Field(alias="posthoc pvalue")
    f_statistics: float = Field(alias="F-statistics")
    posthoc_padj: float = Field(alias="posthoc pvalue adj", ge=0, le=1)
    posthoc_rejected: bool = Field(alias="posthoc rejected", nullable=False)

    posthoc_paired: bool = Field(alias="posthoc Paired", ge=0, le=1)
    posthoc_parametric: bool = Field(alias="posthoc Parametric", ge=0, le=1)
    posthoc_dof: float = Field(alias="posthoc dof", ge=0)
    posthoc_tail: str = Field(alias="posthoc tail")
    posthoc_BF10: float = Field(alias="posthoc BF10", coerce=True)
    posthoc_effsize: float = Field(alias="posthoc effsize")
    efftype: str = Field(alias="efftype")


class AncovaSchema(DataFrameModel):
    """
    Schema for the differential regulation analysis results DataFrame.
    """

    # same as in AnovaSchemaMultiGroup
    posthoc_comparison: str = Field(alias="posthoc comparison", nullable=False)
    t_statistics: float = Field(alias="posthoc T-Statistics")
    posthoc_pvalue: float = Field(alias="posthoc pvalue")
    f_statistics: float = Field(alias="F-statistics")
    posthoc_padj: float = Field(alias="posthoc pvalue adj", ge=0, le=1)
    posthoc_rejected: bool = Field(alias="posthoc rejected", nullable=False)

    coef: float = Field()
    std_err: float = Field(alias="std err")
    conf_int_low: float = Field(alias="Conf. Int. Low")
    conf_int_upp: float = Field(alias="Conf. Int. Upp.")
