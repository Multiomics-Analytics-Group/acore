import pandas as pd


def eta_squared(aov: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the effect size using Eta-squared.

    :param aov: pandas dataframe with anova results from statsmodels.
    :return: Pandas dataframe with additional Eta-squared column.
    """
    aov["eta_sq"] = "NaN"
    aov["eta_sq"] = aov[:-1]["sum_sq"] / sum(aov["sum_sq"])
    return aov


def omega_squared(aov: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the effect size using Omega-squared.

    :param aov: pandas dataframe with anova results from statsmodels.
    :return: Pandas dataframe with additional Omega-squared column.
    """
    mse = aov["sum_sq"][-1] / aov["df"][-1]
    aov["omega_sq"] = "NaN"
    aov["omega_sq"] = (aov[:-1]["sum_sq"] - (aov[:-1]["df"] * mse)) / (
        sum(aov["sum_sq"]) + mse
    )
    return aov
