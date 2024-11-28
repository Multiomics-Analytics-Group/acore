import unittest

import numpy as np
import pandas as pd
import pingouin as pg
from scipy import stats

import acore.correlation_analysis as ca


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
