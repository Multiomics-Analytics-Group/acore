import warnings

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def _get_annotations(data, annotation_cols, group):
    """Select annotation columns (or index levels) with a positional index and
    return the columns which should not be used as features."""
    index_names = [name for name in data.index.names if name is not None]
    missing_annotation_cols = [
        col
        for col in annotation_cols
        if col not in data.columns and col not in index_names
    ]
    if missing_annotation_cols:
        warnings.warn(
            f"Annotation columns not found in data columns or index: "
            f"{missing_annotation_cols}. They are ignored."
        )
    annotation_cols = [
        col for col in annotation_cols if col not in missing_annotation_cols
    ]
    # positional index, as the group labels used as index are not unique
    annotations = pd.DataFrame(
        {
            col: (
                data[col].to_numpy()
                if col in data.columns
                else data.index.get_level_values(col)
            )
            for col in annotation_cols
        },
        index=pd.RangeIndex(len(data)),
    )
    # annotation columns are not used as features
    non_feature_cols = set(annotation_cols).intersection(data.columns)
    return annotations, non_feature_cols.difference([group])


def run_pca(
    data: pd.DataFrame,
    drop_cols: list[str] | None = None,
    group: str = "group",
    annotation_cols: list[str] | None = None,
    components: int = 2,
    dropna: bool = True,
) -> tuple[tuple[pd.DataFrame, pd.DataFrame, np.ndarray | list[float]], dict[str, str]]:
    """
    Performs principal component analysis and returns the values of each component for each sample
    and each protein, and the loadings for each protein.

    For information visit
    https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
                           Default is [] if None.
    :param str group: column label containing group identifiers.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation.
                                 Default is [] if None. Can be columns or index
                                 level names; missing ones are ignored with a warning.
    :param int components: number of components to keep.
    :param bool dropna: if True removes all columns with any missing values.
    :return: tuple: 1) three pandas dataframes: components, loadings and variance; 2)
             xaxis and yaxis titles with components loadings for plotly.

    Example::

        result = run_pca(data, drop_cols=['sample', 'subject'], group='group',
                         components=2, dropna=True)
    """

    if annotation_cols is None:
        annotation_cols = []
    if drop_cols is None:
        drop_cols = []
    np.random.seed(112736)
    var_exp = []
    args = {}
    if data.empty:
        raise ValueError("Dataframe is empty.")

    df = data.copy()
    annotations, non_feature_cols = _get_annotations(data, annotation_cols, group)
    drop_cols_int = list((set(drop_cols) | non_feature_cols).intersection(df.columns))
    if len(drop_cols_int) > 0:
        df = df.drop(drop_cols_int, axis=1)

    y = df[group].tolist()
    df = df.set_index(group)
    df = df.select_dtypes(["number"])
    if dropna:
        df = df.dropna(axis=1)
    X = df.values

    if X.size > 0 and X.shape[1] > components:
        pca = PCA(n_components=components)
        X = pca.fit_transform(X)
        var_exp = pca.explained_variance_ratio_
        loadings = pd.DataFrame(pca.components_.transpose())
        loadings.index = df.columns
        # ? apply calculation before transposition on columns (for performance)
        values = {
            index: np.sqrt(np.power(row, 2).sum()) for index, row in loadings.iterrows()
        }
        loadings["value"] = loadings.index.map(values.get)
        loadings = loadings.sort_values(by="value", ascending=False)
        args = {
            "x": "PC1" + f" ({var_exp[0]:.2f})",
            "y": "PC2" + f" ({var_exp[1]:.2f})",
            "value": "feature_communality",
            "group": "group",
        }
        if components == 2:
            pca_cols = ["x", "y"]
            loadings.columns = ["x", "y", "value"]
        if components > 2:
            args.update({"z": "PC3" + f" ({var_exp[2]:.2f})"})
            pca_cols = ["x", "y", "z"] + [str(i) for i in range(3, components)]
            loading_cols = [str(i) for i in loadings.columns[3:]]
            loadings.columns = ["x", "y", "z"] + loading_cols
        resultDf = pd.DataFrame(X, columns=pca_cols)
        resultDf.insert(0, "group", y)
        resultDf = pd.concat([resultDf, annotations], axis=1)

    return (resultDf, loadings, var_exp), args


def run_tsne(
    data: pd.DataFrame,
    drop_cols: list[str] | None = None,
    group: str = "group",
    annotation_cols: list[str] | None = None,
    components: int = 2,
    perplexity: int = 40,
    max_iter: int = 1000,
    init: str = "pca",
    dropna: bool = True,
    method: str = "barnes_hut",
    random_state: int | None = None,
) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    """
    Performs t-distributed Stochastic Neighbor Embedding analysis.

    For more information visit
    https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
    Default is [] if None.
    :param str group: column label containing group identifiers.
    :param int components: dimension of the embedded space.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation.
                                 Default is [] if None. Can be columns or index
                                 level names; missing ones are ignored with a warning.
    :param int perplexity: related to the number of nearest neighbors that is used
                           in other manifold learning algorithms.
                           Consider selecting a value between 5 and 50.
    :param int max_iter: maximum number of iterations for the optimization (at least 250).
    :param str init: initialization of embedding ('random', 'pca' or
                     numpy array of shape n_samples x n_components).
    :param bool dropna: if True removes all columns with any missing values.
    :param str method: gradient calculation algorithm ('barnes_hut' or 'exact').
                       'barnes_hut' is an approximation computed in float32 whose
                       results can differ between CPU architectures (e.g. arm64
                       vs x64) for small or ill-conditioned data. 'exact' is
                       slower (quadratic) but reproducible across platforms.
    :param int random_state: seed for reproducible results (used with
                             init='random' and the optimization).
    :return: Two dictionaries:
                1) pandas dataframe with embedding vectors,
                2) xaxis and yaxis titles for plotly.

    Example::

        result = run_tsne(data,
                          drop_cols=['sample', 'subject'],
                          group='group',
                          components=2,
                          perplexity=40,
                          max_iter=1000,
                          init='pca',
                          dropna=True
                        )
    """
    if annotation_cols is None:
        annotation_cols = []
    if drop_cols is None:
        drop_cols = []
    result = {}
    args = {}
    df = data.copy()
    annotations, non_feature_cols = _get_annotations(data, annotation_cols, group)
    df = df.drop(
        list((set(drop_cols) | non_feature_cols).intersection(df.columns)), axis=1
    )
    df = df.set_index(group)
    if dropna:
        df = df.dropna(axis=1)
    df = df.select_dtypes(["number"])
    X = df.values
    y = df.index
    if X.size > 0:
        tsne = TSNE(
            n_components=components,
            verbose=0,
            perplexity=perplexity,
            max_iter=max_iter,
            init=init,
            method=method,
            random_state=random_state,
        )
        X = tsne.fit_transform(X)
        args = {"x": "t-SNE1", "y": "t-SNE2"}
        embedding_cols = ["x", "y"]
        if components > 2:
            args.update({"z": "t-SNE3"})
            embedding_cols = ["x", "y", "z"] + [str(i) for i in range(3, components)]
        resultDf = pd.DataFrame(X, columns=embedding_cols)
        resultDf.insert(0, "group", y)
        resultDf = pd.concat([resultDf, annotations], axis=1)
        result["tsne"] = resultDf
    return result, args


def run_umap(
    data: pd.DataFrame,
    drop_cols: list[str] | None = None,
    group: str = "group",
    annotation_cols: list[str] | None = None,
    n_neighbors: int = 10,
    min_dist: float = 0.3,
    metric: str = "cosine",
    dropna: bool = True,
) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    """
    Performs Uniform Manifold Approximation and Projection.

    For more information vist https://umap-learn.readthedocs.io.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
                           Default is [] if None.
    :param str group: column label containing group identifiers.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation.
                                 Default is [] if None. Can be columns or index
                                 level names; missing ones are ignored with a warning.
    :param int n_neighbors: number of neighboring points used
                            in local approximations of manifold structure.
    :param float min_dist: controls how tightly the embedding is allowed compress points together.
    :param str metric: metric used to measure distance in the input space.
    :param bool dropna: if True removes all columns with any missing values.
    :return: Two dictionaries:
                1) pandas dataframe with embedding of the training data in low-dimensional space,
                2) xaxis and yaxis titles for plotly.

    Example::

        result = run_umap(data,
                          drop_cols=['sample', 'subject'],
                          group='group',
                          n_neighbors=10,
                          min_dist=0.3,
                          metric='cosine',
                          dropna=True
                        )
    """
    if annotation_cols is None:
        annotation_cols = []
    if drop_cols is None:
        drop_cols = []
    np.random.seed(1145536)
    result = {}
    args = {}
    df = data.copy()
    annotations, non_feature_cols = _get_annotations(data, annotation_cols, group)
    df = df.drop(
        list((set(drop_cols) | non_feature_cols).intersection(df.columns)), axis=1
    )
    df = df.set_index(group)
    if dropna:
        df = df.dropna(axis=1)
    df = df.select_dtypes(["number"])
    X = df.values
    y = df.index

    if not X.size:
        return result, args

    import umap  # slow import (numba compilation), so only import when needed

    X = umap.UMAP(
        n_neighbors=n_neighbors, min_dist=min_dist, metric=metric
    ).fit_transform(X)
    args = {"x": "UMAP1", "y": "UMAP2"}
    embedding_cols = ["x", "y"] + [str(i) for i in range(2, X.shape[1])]
    resultDf = pd.DataFrame(X, columns=embedding_cols)
    resultDf.insert(0, "group", y)
    resultDf = pd.concat([resultDf, annotations], axis=1)
    result["umap"] = resultDf

    return result, args
