import unittest
from scipy import stats
import acore.enrichment_analysis as ea


class TestRunFisher(unittest.TestCase):
    def test_run_fisher(self):
        group1 = [10, 5]
        group2 = [8, 12]
        alternative = 'two-sided'

        expected_odds, expected_pvalue = stats.fisher_exact([[10, 5], [8, 12]], alternative)

        result = ea.run_fisher(group1, group2, alternative=alternative)

        self.assertEqual(result[0], expected_odds)
        self.assertEqual(result[1], expected_pvalue)


class TestRunKolmogorovSmirnov(unittest.TestCase):
    def test_run_kolmogorov_smirnov(self):
        dist1 = [1, 2, 3, 4, 5]
        dist2 = [1, 2, 3, 4, 6]
        alternative = 'two-sided'

        expected_result = stats.ks_2samp(dist1, dist2, alternative=alternative, mode='auto')

        result = ea.run_kolmogorov_smirnov(dist1, dist2, alternative=alternative)

        self.assertEqual(result[0], expected_result.statistic)
        self.assertEqual(result[1], expected_result.pvalue)


if __name__ == '__main__':
    unittest.main()
