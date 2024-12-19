import unittest

import numpy as np
import pandas as pd
from scipy import stats

import acore.enrichment_analysis as ea


class TestRunFisher(unittest.TestCase):
    def test_run_fisher(self):
        group1 = [10, 5]
        group2 = [8, 12]
        alternative = "two-sided"

        expected_odds, expected_pvalue = stats.fisher_exact(
            [[10, 5], [8, 12]], alternative
        )

        result = ea.run_fisher(group1, group2, alternative=alternative)

        self.assertEqual(result[0], expected_odds)
        self.assertEqual(result[1], expected_pvalue)


class TestRunKolmogorovSmirnov(unittest.TestCase):
    def test_run_kolmogorov_smirnov(self):
        dist1 = [1, 2, 3, 4, 5]
        dist2 = [1, 2, 3, 4, 6]
        alternative = "two-sided"

        expected_result = stats.ks_2samp(
            dist1, dist2, alternative=alternative, mode="auto"
        )

        result = ea.run_kolmogorov_smirnov(dist1, dist2, alternative=alternative)

        self.assertEqual(result[0], expected_result.statistic)
        self.assertEqual(result[1], expected_result.pvalue)


def test__annotate_features():
    expected = pd.Series(
        [
            "foreground",
            "foreground",
            "background",
            "foreground",
            "background",
            "background",
            np.nan,
        ]
    )

    features = pd.Series(["G1", "G2", "G3", "G4", "G5", "G6", "G9"])
    in_foreground = ["G1", "G2", "G4"]
    in_background = ["G3", "G5", "G6"]
    actual = ea._annotate_features(features, in_foreground, in_background)
    pd.testing.assert_series_equal(expected, actual)


def test_run_regulation_enrichment():
    """Integration test for run_regulation_enrichment. Indirectly tests
    run_enrichment from enrichment_analysis module."""
    annotation = {
        "annotation": ["path1", "path1", "path1", "path2", "path2", "path3", "path3"],
        "identifier": ["gene1", "gene2", "gene3", "gene1", "gene5", "gene6", "gene9"],
        "source": ["GO", "GO", "GO", "GO_P", "GO_P", "GO_P", "GO_P"],
    }
    annotation = pd.DataFrame(annotation)
    regulation_res = {
        "identifier": ["gene1", "gene2", "gene3", "gene4", "gene5", "gene6"],
        "rejected": [True, True, False, False, True, True],
    }
    regulation_res = pd.DataFrame(regulation_res)

    actual = ea.run_regulation_enrichment(
        regulation_data=regulation_res,
        annotation=annotation,
        min_detected_in_set=1,
    )

    expected = pd.DataFrame(
        {
            "terms": ["path1", "path2", "path3"],
            "identifiers": ["gene1,gene2", "gene1,gene5", "gene6"],
            "foreground": [2, 2, 1],
            "background": [1, 0, 0],
            "foreground_pop": [4, 4, 4],
            "background_pop": [6, 6, 6],
            "pvalue": [1.0, 0.4666666666666667, 1.0],
            "padj": [1.0, 1.0, 1.0],
            "rejected": [False, False, False],
        }
    )
    assert expected.equals(actual)


if __name__ == "__main__":
    unittest.main()
