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


if __name__ == "__main__":
    unittest.main()
