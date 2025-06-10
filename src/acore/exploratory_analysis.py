import numpy as np
import pandas as pd
import scipy.stats
import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def calculate_coefficient_variation(values: np.ndarray) -> np.ndarray:
    """
    Compute the coefficient of variation, the ratio of the biased standard
    deviation to the mean, in percentage. For more information
    visit https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.variation.html.

    :param numpy.ndarray values: numpy array of log2 transformed values
    :return: The calculated variation along rows.
    :rtype: numpy.ndarray

    Example::

        result = calculate_coefficient_variation()
    """
    cv = scipy.stats.variation(values.apply(lambda x: np.power(2, x)).values) * 100

    return cv


def get_coefficient_variation(data, drop_columns, group, columns=["name", "y"]):
    """
    Extracts the coefficients of variation in each group.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_columns: column labels to be dropped from the dataframe
    :param str group: column label containing group identifiers.
    :param list columns: names to use for the variable column(s), and for the value column(s)
    :return: Pandas dataframe with columns 'name' (protein identifier),
             'x' (coefficient of variation), 'y' (mean) and 'group'.

    Example::

        result = get_coefficient_variation(data, drop_columns=['sample', 'subject'], group='group')
    """
    df = data.copy()
    formated_df = df.drop(drop_columns, axis=1)
    cvs = formated_df.groupby(group).apply(func=calculate_coefficient_variation)
    cols = formated_df.set_index(group).columns.tolist()
    cvs_df = pd.DataFrame()
    for i in cvs.index:
        gcvs = cvs[i].tolist()
        ints = formated_df.set_index(group).mean().values.tolist()
        tdf = pd.DataFrame(data={"name": cols, "x": gcvs, "y": ints})
        tdf[group] = i

        if cvs_df.empty:
            cvs_df = tdf.copy()
        else:
            cvs_df = pd.concat([cvs_df, tdf])

    return cvs_df


def extract_number_missing(data, min_valid, drop_cols=["sample"], group="group"):
    """
    Counts how many valid values exist in each column and filters column labels with more
    valid values than the minimum threshold defined.

    :param data: pandas DataFrame with group as rows and protein identifier as column.
    :param str group: column label containing group identifiers.
                      If None, number of valid values is counted across all samples,
                      otherwise is counted per unique group identifier.
    :param int min_valid: minimum number of valid values to be filtered.
    :param list drop_columns: column labels to be dropped.
    :return: List of column labels above the threshold.

    Example::

        result = extract_number_missing(data, min_valid=3, drop_cols=['sample'], group='group')
    """
    if group is None:
        groups = data.loc[:, data.notnull().sum(axis=0) >= min_valid]
    else:
        groups = data.copy()
        if len(set(drop_cols).intersection(groups.columns.tolist())) == len(drop_cols):
            groups = groups.drop(drop_cols, axis=1)
        groups = groups.set_index(group).notnull().groupby(level=0).sum()
        groups = groups[groups >= min_valid]

    groups = groups.dropna(how="all", axis=1)
    return groups.columns.unique().tolist()


def extract_percentage_missing(
    data, missing_max, drop_cols=["sample"], group="group", how="all"
):
    """
    Extracts ratio of missing/valid values in each column and filters column labels with
    lower ratio than the minimum threshold defined.

    :param data: pandas dataframe with group as rows and protein identifier as column.
    :param str group: column label containing group identifiers.
                      If None, ratio is calculated across all samples,
                      otherwise is calculated per unique group identifier.
    :param float missing_max: maximum ratio of missing/valid values to be filtered.
    :param str how: define if labels with a higher percentage of missing values than the threshold
                    in any group ('any') or in all groups ('all') should be filtered
    :return: List of column labels below the threshold.

    Example::
        result = extract_percentage_missing(data, missing_max=0.3,
                                            drop_cols=['sample'], group='group')
    """
    if group is None:
        groups = data.loc[:, data.isnull().mean() <= missing_max].columns
    else:
        groups = data.copy()
        groups = groups.drop(drop_cols, axis=1)
        groups = groups.set_index(group)
        groups = groups.isnull().groupby(level=0).mean()
        groups = groups[groups <= missing_max]
        if how == "all":
            groups = groups.dropna(how="all", axis=1).columns.unique().tolist()
        elif how == "any":
            groups = groups.dropna(how="any", axis=1).columns.unique().tolist()
        else:
            if how in groups.index:
                groups = groups.loc[how, :].dropna().index.unique().tolist()

    return groups


def run_pca(
    data,
    drop_cols=["sample", "subject"],
    group="group",
    annotation_cols=["sample"],
    components=2,
    dropna=True,
):
    """
    Performs principal component analysis and returns the values of each component for each sample
    and each protein, and the loadings for each protein.

    For information visit
    https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
    :param str group: column label containing group identifiers.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation
    :param int components: number of components to keep.
    :param bool dropna: if True removes all columns with any missing values.
    :return: tuple: 1) three pandas dataframes: components, loadings and variance; 2)
             xaxis and yaxis titles with components loadings for plotly.

    Example::

        result = run_pca(data, drop_cols=['sample', 'subject'], group='group',
                         components=2, dropna=True)
    """

    np.random.seed(112736)
    var_exp = []
    args = {}
    if data.empty:
        raise ValueError("Dataframe is empty.")

    df = data.copy()
    annotations = pd.DataFrame()
    if annotation_cols is not None:
        if len(list(set(annotation_cols).intersection(data.columns))) > 0:
            annotations = data.set_index(group)[annotation_cols]
    drop_cols_int = list(set(drop_cols).intersection(df.columns))
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
        values = {
            index: np.sqrt(np.power(row, 2).sum()) for index, row in loadings.iterrows()
        }
        loadings["value"] = loadings.index.map(values.get)
        loadings = loadings.sort_values(by="value", ascending=False)
        args = {
            "x_title": "PC1" + " ({0:.2f})".format(var_exp[0]),
            "y_title": "PC2" + " ({0:.2f})".format(var_exp[1]),
            "group": "group",
        }
        if components == 2:
            resultDf = pd.DataFrame(X, index=y, columns=["x", "y"])
            resultDf = resultDf.assign(**annotations)
            resultDf = resultDf.reset_index()
            resultDf.columns = ["group", "x", "y"] + annotation_cols

            loadings.columns = ["x", "y", "value"]
        if components > 2:
            args.update({"z_title": "PC3" + " ({0:.2f})".format(var_exp[2])})
            resultDf = pd.DataFrame(X, index=y)
            resultDf = resultDf.assign(**annotations)
            resultDf = resultDf.reset_index()
            pca_cols = []
            loading_cols = []
            if components > 3:
                pca_cols = [str(i) for i in resultDf.columns[4:]]
                loading_cols = [str(i) for i in loadings.columns[3:]]

            resultDf.columns = ["group", "x", "y", "z"] + pca_cols
            loadings.columns = ["x", "y", "z"] + loading_cols

    return (resultDf, loadings, var_exp), args


def run_tsne(
    data,
    drop_cols=["sample", "subject"],
    group="group",
    annotation_cols=["sample"],
    components=2,
    perplexity=40,
    max_iter=1000,
    init="pca",
    dropna=True,
):
    """
    Performs t-distributed Stochastic Neighbor Embedding analysis.

    For more information visit
    https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
    :param str group: column label containing group identifiers.
    :param int components: dimension of the embedded space.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation
    :param int perplexity: related to the number of nearest neighbors that is used
                           in other manifold learning algorithms.
                           Consider selecting a value between 5 and 50.
    :param int max_iter: maximum number of iterations for the optimization (at least 250).
    :param str init: initialization of embedding ('random', 'pca' or
                     numpy array of shape n_samples x n_components).
    :param bool dropna: if True removes all columns with any missing values.
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
    result = {}
    args = {}
    df = data.copy()
    if len(set(drop_cols).intersection(df.columns)) == len(drop_cols):
        df = df.drop(drop_cols, axis=1)
    df = df.set_index(group)
    if dropna:
        df = df.dropna(axis=1)
    df = df.select_dtypes(["number"])
    X = df.values
    y = df.index
    annotations = pd.DataFrame()
    if annotation_cols is not None:
        if len(list(set(annotation_cols).intersection(data.columns))) > 0:
            annotations = data[annotation_cols]
    if X.size > 0:
        tsne = TSNE(
            n_components=components,
            verbose=0,
            perplexity=perplexity,
            max_iter=max_iter,
            init=init,
        )
        X = tsne.fit_transform(X)
        args = {"x_title": "C1", "y_title": "C2"}
        if components == 2:
            resultDf = pd.DataFrame(X, index=y, columns=["x", "y"])
            resultDf = resultDf.reset_index()
            resultDf.columns = ["group", "x", "y"]
        if components > 2:
            args.update({"z_title": "C3"})
            resultDf = pd.DataFrame(X, index=y)
            resultDf = resultDf.reset_index()
            cols = []
            if len(components) > 4:
                cols = resultDf.columns[4:]
            resultDf.columns = ["group", "x", "y", "z"] + cols
        resultDf = resultDf.join(annotations)
        result["tsne"] = resultDf
    return result, args


def run_umap(
    data,
    drop_cols=["sample", "subject"],
    group="group",
    annotation_cols=["sample"],
    n_neighbors=10,
    min_dist=0.3,
    metric="cosine",
    dropna=True,
):
    """
    Performs Uniform Manifold Approximation and Projection.

    For more information vist https://umap-learn.readthedocs.io.

    :param data: pandas dataframe with samples as rows and protein identifiers as columns
                 (with additional columns 'group', 'sample' and 'subject').
    :param list drop_cols: column labels to be dropped from the dataframe.
    :param str group: column label containing group identifiers.
    :param list annotation_cols: list of columns to be added in the scatter plot annotation
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
    np.random.seed(1145536)
    result = {}
    args = {}
    df = data.copy()
    if len(set(drop_cols).intersection(df.columns)) == len(drop_cols):
        df = df.drop(drop_cols, axis=1)
    df = df.set_index(group)
    if dropna:
        df = df.dropna(axis=1)
    df = df.select_dtypes(["number"])
    X = df.values
    y = df.index

    annotations = pd.DataFrame()
    if annotation_cols is not None:
        if len(list(set(annotation_cols).intersection(data.columns))) > 0:
            annotations = data[annotation_cols]

    if X.size:
        X = umap.UMAP(
            n_neighbors=n_neighbors, min_dist=min_dist, metric=metric
        ).fit_transform(X)
        args = {"x_title": "C1", "y_title": "C2"}
        resultDf = pd.DataFrame(X, index=y)
        resultDf = resultDf.reset_index()
        cols = []
        if len(resultDf.columns) > 3:
            cols = resultDf.columns[3:]
        resultDf.columns = ["group", "x", "y"] + cols
        resultDf = resultDf.join(annotations)
        result["umap"] = resultDf

    return result, args
