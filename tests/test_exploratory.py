import unittest
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
import acore.exploratory_analysis as ea

class TestExtractMissing(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame({
            'group': ['A', 'A', 'B', 'B'],
            'protein1': [1, 2, None, 4],
            'protein2': [5, None, 7, 8],
            'protein3': [9, 10, 11, 12]
        })

    def test_extract_number_missing(self):
        min_valid = 2
        expected_result = ['protein1', 'protein3']

        result = ea.extract_number_missing(self.data, min_valid, drop_cols=['sample'], group='group')

        self.assertEqual(result, expected_result)

    def test_extract_percentage_missing(self):
        missing_max = 0.3
        expected_result = ['protein2']

        result = ea.extract_percentage_missing(self.data, missing_max, drop_cols=['sample'], group='group')

        self.assertEqual(result, expected_result)



class TestDimensionalityReduction(unittest.TestCase):
    def setUp(self):
        self.data = pd.DataFrame({
            'group': ['A', 'A', 'B', 'B'],
            'protein1': [1, 2, None, 4],
            'protein2': [5, None, 7, 8],
            'protein3': [9, 10, 11, 12]
        })

    def test_run_pca(self):
        components = 2
        expected_result = (
            pd.DataFrame({
                'group': ['A', 'A', 'B', 'B'],
                'x': [0.212769, 0.219452, -0.521191, -0.214031],
                'y': [-0.116136, -0.116432, 0.276313, 0.113255],
                'protein1': [1, 2, None, 4],
                'protein2': [5, None, 7, 8],
                'protein3': [9, 10, 11, 12]
            }),
            pd.DataFrame({
                'x': [0.58497, 0.655434, 0.47958],
                'y': [0.395222, 0.313202, 0.863743],
                'value': [0.851531, 0.781998, 0.18541]
            }),
            [0.548550, 0.337207, 0.114243]
        )

        result, args = run_pca(self.data, drop_cols=['sample', 'subject'], group='group', components=components, dropna=True)

        self.assertEqual(result, expected_result[0])
        self.assertEqual(args, {'x_title': 'PC1 (0.55)', 'y_title': 'PC2 (0.34)', 'group': 'group'})

    def test_run_tsne(self):
        components = 2
        expected_result = (
            {
                'tsne': pd.DataFrame({
                    'group': ['A', 'A', 'B', 'B'],
                    'x': [0.212769, 0.219452, -0.521191, -0.214031],
                    'y': [-0.116136, -0.116432, 0.276313, 0.113255]
                })
            },
            {'x_title': 'C1', 'y_title': 'C2'}
        )

        result, args = run_tsne(self.data, drop_cols=['sample', 'subject'], group='group', components=components, perplexity=40, n_iter=1000, init='pca', dropna=True)

        self.assertEqual(result, expected_result[0])
        self.assertEqual(args, expected_result[1])

    def test_run_umap(self):
        expected_result = (
            {
                'umap': pd.DataFrame({
                    'group': ['A', 'A', 'B', 'B'],
                    'x': [0.234967, 0.239225, -0.562752, -0.231992],
                    'y': [-0.121199, -0.121862, 0.284414, 0.117647]
                })
            },
            {'x_title': 'C1', 'y_title': 'C2'}
        )

        result, args = run_umap(self.data, drop_cols=['sample', 'subject'], group='group', n_neighbors=10, min_dist=0.3, metric='cosine', dropna=True)

        self.assertEqual(result, expected_result[0])
        self.assertEqual(args, expected_result[1])

if __name__ == '__main__':
    unittest.main()
