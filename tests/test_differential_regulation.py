import unittest

import pandas as pd
import pingouin as pg

import acore.differential_regulation as dr


class TestCalculateTtest(unittest.TestCase):
    def setUp(self):
        self.data = {
            "subject": [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
            "sample": [
                "S1",
                "S2",
                "S3",
                "S4",
                "S5",
                "S6",
                "S7",
                "S8",
                "S9",
                "S10",
                "S11",
                "S12",
            ],
            "group": ["A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"],
            "protein": [
                1.4,
                6.2,
                9.03,
                2.3,
                7.4,
                14.01,
                4.5,
                9.6,
                10.4,
                7.5,
                11.6,
                16.4,
            ],
            "group2": ["a", "a", "a", "a", "a", "a", "b", "b", "b", "b", "b", "b"],
        }
        self.df = pd.DataFrame(self.data)

    def test_ttest_ind(self):
        data = {
            "group1": [1.3, 1.2, 0.2, 2.2, 3.3, 8.43],
            "group2": [1.4, 6.2, 7.3, 0.4, 1.5, 0.6],
        }
        df = pd.DataFrame(data)

        condition1 = "group1"
        condition2 = "group2"
        is_logged = True
        non_par = False

        expected_result = pg.ttest(df[condition1], df[condition2], paired=False)
        result = dr.calculate_ttest(
            df, condition1, condition2, is_logged=is_logged, non_par=non_par
        )

        self.assertAlmostEqual(result[0], expected_result["T"].values[0])
        self.assertAlmostEqual(result[1], expected_result["p-val"].values[0])

    def test_mann_whitney(self):
        data = {"group1": [1, 2, 3, 4, 5], "group2": [2, 4, 6, 8, 10]}
        df = pd.DataFrame(data)

        condition1 = "group1"
        condition2 = "group2"
        is_logged = True
        non_par = True

        expected_result = pg.mwu(df[condition1], df[condition2])
        result = dr.calculate_ttest(
            df, condition1, condition2, is_logged=is_logged, non_par=non_par
        )

        self.assertAlmostEqual(result[0], expected_result["U-val"].values[0])
        self.assertAlmostEqual(result[1], expected_result["p-val"].values[0])

    def test_calculate_anova(self):
        column = "protein"
        group = "group"

        expected_result = pg.anova(data=self.df, dv=column, between=group)
        expected_t, expected_df1, expected_df2, expected_pvalue = expected_result[
            ["F", "ddof1", "ddof2", "p-unc"]
        ].values[0]

        result = dr.calculate_anova(self.df, column, group=group)

        self.assertEqual(result[1], expected_df1)
        self.assertEqual(result[2], expected_df2)
        self.assertEqual(result[3], expected_t)
        self.assertEqual(result[4], expected_pvalue)

    def test_calculate_ancova(self):
        column = "protein"
        group = "group"
        covariates = []

        expected_result = pg.ancova(
            data=self.df, dv=column, between=group, covar=covariates
        )
        expected_t, expected_df, expected_pvalue = expected_result.loc[
            expected_result["Source"] == group, ["F", "DF", "p-unc"]
        ].values[0]

        result = dr.calculate_ancova(
            self.df, column, group=group, covariates=covariates
        )

        self.assertEqual(result[1], expected_df)
        self.assertEqual(result[2], expected_df)
        self.assertEqual(result[3], expected_t)
        self.assertEqual(result[4], expected_pvalue)

    def test_calculate_repeated_measures_anova(self):
        """Source     SS  DF      MS         F     p-unc       ng2  eps
        0  group   1.50   1   1.500  0.087642  0.795106  0.032609  1.0
        1  Error  34.23   2  17.115       NaN       NaN       NaN  NaN"""
        column = "protein"
        subject = "subject"
        within = "group"

        expected_result = pg.rm_anova(
            data=self.df,
            dv=column,
            within=within,
            subject=subject,
            detailed=True,
            correction=True,
        )
        expected_t, expected_pvalue = expected_result.loc[
            0, ["F", "p-unc"]
        ].values.tolist()
        expected_df1, expected_df2 = expected_result["DF"]
        result = dr.calculate_repeated_measures_anova(
            self.df, column, subject=subject, within=within
        )

        self.assertEqual(result[1], expected_df1)
        self.assertEqual(result[2], expected_df2)
        self.assertEqual(result[3], expected_t)
        self.assertEqual(result[4], expected_pvalue)

    def test_calculate_mixed_anova(self):
        column = "protein"
        subject = "subject"
        within = "group"
        between = "group2"

        expected_result = pg.mixed_anova(
            data=self.df,
            dv=column,
            within=within,
            between=between,
            subject=subject,
            correction=True,
        )
        expected_result["identifier"] = column
        expected_result = expected_result[
            ["identifier", "DF1", "DF2", "F", "p-unc", "Source"]
        ]

        result = dr.calculate_mixed_anova(
            self.df, column, subject=subject, within=within, between=between
        )

        self.assertEqual(result.values.tolist(), expected_result.values.tolist())


class TestDropColsOptional(unittest.TestCase):
    """Tests that drop_cols=None and drop_cols=[] are handled correctly in all functions."""

    def setUp(self):
        # Two-group dataset without extra non-protein columns
        # (demonstrating that drop_cols=None works when the df is already clean)
        self.data_2groups = {
            "subject": [1, 2, 3, 4, 1, 2, 3, 4],
            "group": ["A", "A", "A", "A", "B", "B", "B", "B"],
            "protein": [1.4, 2.3, 4.5, 7.5, 6.2, 7.4, 9.6, 11.6],
        }
        self.df_2groups = pd.DataFrame(self.data_2groups)

        # Three-group dataset without extra non-protein columns
        self.data_3groups = {
            "subject": [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4],
            "group": ["A", "A", "A", "A", "B", "B", "B", "B", "C", "C", "C", "C"],
            "protein": [1.4, 2.3, 4.5, 7.5, 6.2, 7.4, 9.6, 11.6, 9.0, 14.0, 10.4, 16.4],
        }
        self.df_3groups = pd.DataFrame(self.data_3groups)

    def test_run_ttest_drop_cols_none(self):
        """run_ttest should not crash when drop_cols=None."""
        result = dr.run_ttest(
            self.df_2groups,
            condition1="A",
            condition2="B",
            drop_cols=None,
            group="group",
            subject="subject",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_ttest_drop_cols_empty(self):
        """run_ttest should not crash when drop_cols=[]."""
        result = dr.run_ttest(
            self.df_2groups,
            condition1="A",
            condition2="B",
            drop_cols=[],
            group="group",
            subject="subject",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_ttest_drop_cols_none_matches_empty(self):
        """run_ttest results with drop_cols=None and drop_cols=[] should be identical."""
        result_none = dr.run_ttest(
            self.df_2groups,
            condition1="A",
            condition2="B",
            drop_cols=None,
            group="group",
            subject="subject",
        )
        result_empty = dr.run_ttest(
            self.df_2groups,
            condition1="A",
            condition2="B",
            drop_cols=[],
            group="group",
            subject="subject",
        )
        pd.testing.assert_frame_equal(result_none, result_empty)

    def test_run_anova_two_groups_drop_cols_none(self):
        """run_anova should not crash when drop_cols=None (2-group case)."""
        result = dr.run_anova(
            self.df_2groups,
            drop_cols=None,
            subject="subject",
            group="group",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_anova_two_groups_drop_cols_empty(self):
        """run_anova should not crash when drop_cols=[] (2-group case)."""
        result = dr.run_anova(
            self.df_2groups,
            drop_cols=[],
            subject="subject",
            group="group",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_anova_three_groups_drop_cols_none(self):
        """run_anova should not crash when drop_cols=None (>2-group case, unpaired)."""
        result = dr.run_anova(
            self.df_3groups,
            drop_cols=None,
            subject=None,
            group="group",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_anova_three_groups_drop_cols_empty(self):
        """run_anova should not crash when drop_cols=[] (>2-group case, unpaired)."""
        result = dr.run_anova(
            self.df_3groups,
            drop_cols=[],
            subject=None,
            group="group",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_anova_paired_three_groups_drop_cols_none(self):
        """run_anova should not crash when drop_cols=None (>2-group paired case).

        This also verifies the fix for the inconsistency where subject was passed
        in drop_cols to run_repeated_measurements_anova (which still needs subject).
        """
        result = dr.run_anova(
            self.df_3groups,
            drop_cols=None,
            subject="subject",
            group="group",
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)

    def test_run_repeated_measurements_anova_drop_cols_none(self):
        """run_repeated_measurements_anova should not crash when drop_cols=None."""
        result = dr.run_repeated_measurements_anova(
            self.df_3groups,
            drop_cols=None,
            subject="subject",
            within="group",
            permutations=0,
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn("identifier", result.columns)


if __name__ == "__main__":
    unittest.main()
