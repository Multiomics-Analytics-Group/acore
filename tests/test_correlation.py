import unittest

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

import acore.correlation_analysis as ca

data = {
    "Uracil": {
        "rapZE227Stop-1": 26.270566729949532,
        "rapZE227Stop-2": 27.293323104311497,
        "rapZE227Stop-3": 26.78811750021429,
        "rapZE227Stop-4": 27.014167948711147,
        "rapZE227Stop-5": 26.363653630201306,
        "rapZE227Stop-6": 26.571155418705217,
        "WT-1": 26.47152940077166,
        "WT-2": 26.32477219688337,
        "WT-3": 26.274162656291463,
        "WT-4": 26.1871201052153,
        "WT-5": 26.24356100844659,
        "WT-6": 26.313000315279105,
    },
    "Adenine": {
        "rapZE227Stop-1": 21.873918881276147,
        "rapZE227Stop-2": 21.982539026967927,
        "rapZE227Stop-3": 21.83259951503248,
        "rapZE227Stop-4": 22.1197774615395,
        "rapZE227Stop-5": 21.95663194590021,
        "rapZE227Stop-6": 21.884808226307076,
        "WT-1": 21.825610422927305,
        "WT-2": 21.6519413723448,
        "WT-3": 21.65132912836582,
        "WT-4": 21.78047518052626,
        "WT-5": 21.701526015895677,
        "WT-6": 21.76966959105173,
    },
    "Hypoxanthine": {
        "rapZE227Stop-1": 23.20057740721078,
        "rapZE227Stop-2": 23.357469437159953,
        "rapZE227Stop-3": 23.16674222421689,
        "rapZE227Stop-4": 23.51672102681317,
        "rapZE227Stop-5": 23.24567700614833,
        "rapZE227Stop-6": 23.121489367219876,
        "WT-1": 22.980021484616568,
        "WT-2": 22.925293532285732,
        "WT-3": 22.885493913531445,
        "WT-4": 22.93842585475622,
        "WT-5": 22.944516999920193,
        "WT-6": 22.850814798397437,
    },
}


class TestCalculateCorrelations(unittest.TestCase):
    def test_pearson_correlation(self):
        x = np.array([1.5, 0.2, 3.3, 4.34, 5.03])
        y = np.array([2.04, 4.9, 3.6, 0.8, 10.9])
        coefficient, pvalue = ca.calculate_correlations(x, y, method="pearson")
        expected_coefficient, expected_pvalue = stats.pearsonr(x, y)
        self.assertAlmostEqual(coefficient, expected_coefficient)
        self.assertAlmostEqual(pvalue, expected_pvalue)

    def test_spearman_correlation(self):
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        coefficient, pvalue = ca.calculate_correlations(x, y, method="spearman")
        expected_coefficient, expected_pvalue = stats.spearmanr(x, y)
        self.assertAlmostEqual(coefficient, expected_coefficient)
        self.assertAlmostEqual(pvalue, expected_pvalue)

    def test_calculate_rm_correlation(self):
        # Sample test data
        df = pg.read_dataset("rm_corr")
        x = "pH"
        y = "PacO2"
        subject = "Subject"
        # Expected output
        expected_result = pg.rm_corr(data=df, x=x, y=y, subject=subject)

        # Call the function
        result = ca.calculate_rm_correlation(df, x, y, subject)

        # Compare the results
        self.assertAlmostEqual(result[2], expected_result["r"].values[0])
        self.assertAlmostEqual(result[3], expected_result["pval"].values[0])
        self.assertEqual(result[4], expected_result["dof"].values[0])


def test_run_efficient_correlation_pearson():
    df = pd.DataFrame(data)
    corr, p = ca.run_efficient_correlation(df, method="pearson")
    for i in range(df.shape[1]):
        for j in range(max(i, df.shape[1]), df.shape[1]):
            print(i, j)
            corr_, p_ = stats.pearsonr(df.iloc[:, i], df.iloc[:, j])
            np.testing.assert_almost_equal(corr[j, i], corr_)
            if i != j:
                # ToDo: Diagonal p-values should be zero, not one
                np.testing.assert_almost_equal(p[j, i], p_)


def test_corr_lower_triangle():
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 7], "C": [6, 8, 9]})
    expected_result = pd.DataFrame(
        {
            "A": {"A": np.nan, "B": 0.9819805060619659, "C": 0.9819805060619656},
            "B": {"A": np.nan, "B": np.nan, "C": 0.9285714285714283},
            "C": {"A": np.nan, "B": np.nan, "C": np.nan},
        }
    )
    result = ca.corr_lower_triangle(df)
    pd.testing.assert_frame_equal(result, expected_result)


if __name__ == "__main__":
    unittest.main()
