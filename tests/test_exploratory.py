import platform
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

        result = ea.extract_number_missing(
            self.data, min_valid, drop_cols=["sample"], group="group"
        )

        self.assertEqual(result, expected_result)

    def test_extract_percentage_missing(self):
        missing_max = 0.2
        expected_result = ["protein1", "protein3"]

        result = ea.extract_percentage_missing(
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
        exp_dict = {"x_title": "PC1 (0.72)", "y_title": "PC2 (0.26)", "group": "group"}
        result, annotation = ea.run_pca(
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

    def test_run_tsne(self):
        components = 2
        # ! MAC ARM64 yields these results (using the similar installation)
        expected_result = pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "x": [-113.341728, 36.020717, 57.2992514, -134.729141],
                "y": [131.169082, -59.616630, 110.601203, -39.086254],
            }
        )
        # ! Ubuntu x64 yields these results:
        if not platform.machine() == "arm64":
            expected_result = pd.DataFrame(
                {
                    "group": ["A", "A", "B", "B"],
                    "x": [
                        -27.375560760498047,
                        -7.526909828186035,
                        -28.731882095336914,
                        -8.872751235961914,
                    ],
                    "y": [
                        1.547648310661316,
                        0.17019487917423248,
                        -18.321491241455078,
                        -19.713821411132812,
                    ],
                }
            )
        result, _ = ea.run_tsne(
            self.data,
            drop_cols=["sample", "subject"],
            group="group",
            components=components,
            perplexity=3,
            max_iter=1000,
            init="pca",
            dropna=True,
        )
        result = result["tsne"]
        pd.testing.assert_frame_equal(
            result, expected_result, check_exact=False, check_dtype=False
        )

    def test_run_umap(self):
        _ = ea.run_umap(
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
