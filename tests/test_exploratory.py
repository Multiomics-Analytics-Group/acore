import unittest

import pandas as pd

import acore.exploratory_analysis as ea


class TestExtractMissing(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "sample": ["S1", "S2", "S3", "S4"],
                "protein1": [1.4, 2.2, None, 4.2],
                "protein2": [5.6, None, None, 8.1],
                "protein3": [9.1, 10.01, 11.2, 12.9],
            }
        )

    def test_extract_number_missing(self):
        min_valid = 2
        expected_result = ["protein1", "protein3"]

        result = ea.missing.extract_number_missing(
            self.data, min_valid, drop_cols=["sample"], group="group"
        )

        self.assertEqual(result, expected_result)

    def test_extract_percentage_missing(self):
        missing_max = 0.2
        expected_result = ["protein1", "protein3"]

        result = ea.missing.extract_percentage_missing(
            self.data, missing_max, drop_cols=["sample"], group="group"
        )

        self.assertEqual(result, expected_result)


class TestDimensionalityReduction(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "protein1": [1.4, 2.2, 5.3, 4.2],
                "protein2": [5.6, 0.3, 2.1, 8.1],
                "protein3": [9.1, 10.01, 11.2, 12.9],
            }
        )

    def test_run_pca(self):
        components = 2
        expected_result = (
            pd.DataFrame(
                {
                    "group": ["A", "A", "B", "B"],
                    "x": [0.877154, -3.883458, -1.550791, 4.557095],
                    "y": [-2.815737, -0.420979, 2.278209, 0.958507],
                }
            ),
            pd.DataFrame(
                {
                    "x": [0.958489, 0.092382, 0.269749],
                    "y": [-0.235412, 0.790179, 0.565861],
                    "value": [0.986975, 0.795561, 0.626868],
                },
                index=["protein2", "protein1", "protein3"],
            ),
            pd.Series([0.71760534, 0.26139782]),
        )
        exp_dict = {
            "x": "PC1 (0.72)",
            "y": "PC2 (0.26)",
            "value": "feature_communality",
            "group": "group",
        }
        result, annotation = ea.decomposition.run_pca(
            self.data,
            drop_cols=[],
            annotation_cols=[],
            group="group",
            components=components,
            dropna=True,
        )
        assert annotation == exp_dict
        # results: tuple of 3 DFs: components, loadings, explained_variance
        pd.testing.assert_frame_equal(result[0], expected_result[0], check_exact=False)
        pd.testing.assert_frame_equal(result[1], expected_result[1], check_exact=False)
        pd.testing.assert_series_equal(
            pd.Series(result[2]),
            expected_result[2],
            check_exact=False,
            check_dtype=False,
        )

    def test_run_pca_warns_for_missing_annotation_columns(self):
        with self.assertWarnsRegex(UserWarning, "Annotation columns not found.*sample"):
            (result, _, _), _ = ea.decomposition.run_pca(
                self.data,
                annotation_cols=["sample"],
                group="group",
                components=2,
            )
        self.assertListEqual(result.columns.tolist(), ["group", "x", "y"])

    def test_run_pca_annotation_from_index(self):
        data = self.data.set_index(pd.Index(["S1", "S2", "S3", "S4"], name="Sample ID"))
        (result, _, _), _ = ea.decomposition.run_pca(
            data, annotation_cols=["Sample ID"], group="group", components=2
        )
        self.assertListEqual(result["Sample ID"].tolist(), ["S1", "S2", "S3", "S4"])

    def test_run_tsne(self):
        components = 2
        # method="exact" is used as the default Barnes-Hut approximation (float32)
        # gives different results on arm64 and x64 for this tiny dataset.
        expected_result = pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "x": [40.735077, 207.908371, 23.888744, -143.473709],
                "y": [192.089203, 8.426460, -158.064987, 25.697359],
            }
        )
        # not that

        result, _ = ea.decomposition.run_tsne(
            self.data,
            drop_cols=["sample", "subject"],
            group="group",
            components=components,
            perplexity=3,
            max_iter=1000,
            init="pca",
            dropna=True,
            method="exact",
            random_state=0,
        )
        result = result["tsne"]
        pd.testing.assert_frame_equal(
            result, expected_result, check_exact=False, check_dtype=False
        )

    def test_run_umap(self):
        _ = ea.decomposition.run_umap(
            self.data,
            drop_cols=["sample", "subject"],
            group="group",
            n_neighbors=10,
            min_dist=0.3,
            metric="cosine",
            dropna=True,
        )

        self.assertEqual(1, 1)


if __name__ == "__main__":
    unittest.main()
