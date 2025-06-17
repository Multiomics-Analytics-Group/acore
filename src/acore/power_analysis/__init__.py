import itertools

import numpy as np
import pandas as pd
from statsmodels.stats.power import FTestAnovaPower

from acore import utils


def power_analysis(
    data,
    group="group",
    groups=None,
    alpha=0.05,
    power=0.8,
    dep_var="nobs",
    figure=False,
):
    quantiles = ["25% qtl es", "mean es", "50% qtl es", "75% qtl es"]
    if groups is None:
        groups = data[group].unique().tolist()
    k_groups = len(groups)
    effect_sizes = set()
    for col in data.drop([group], axis=1).columns:
        for g1, g2 in itertools.combinations(groups, 2):
            sample1 = data.loc[data[group] == g1, col].values
            sample2 = data.loc[data[group] == g2, col].values
            eff_size = np.abs(utils.cohens_d(sample1, sample2, ddof=1))
            effect_sizes.add(eff_size)

    summary_eff = []
    if len(effect_sizes):
        effect_sizes = list(effect_sizes)
        summary_eff = [
            np.percentile(effect_sizes, 25),
            np.mean(effect_sizes),
            np.percentile(effect_sizes, 50),
            np.percentile(effect_sizes, 75),
        ]

    analysis = FTestAnovaPower()
    sample_sizes = np.array(range(3, 150))
    power_list = []
    labels = []
    samples = []
    for ii, es in enumerate(summary_eff):
        p = analysis.power(es, sample_sizes, alpha, k_groups)
        labels.extend(["%s = %4.2F" % (quantiles[ii], es)] * len(p))
        power_list.extend(p)
        samples.extend(sample_sizes)

    power_df = pd.DataFrame(
        data=list(zip(power_list, samples, labels)),
        columns=["power", "#samples", "labels"],
    )
    sample_size = analysis.solve_power(
        summary_eff[1], power=power, alpha=alpha, k_groups=k_groups
    )

    return (sample_size, power_df)
