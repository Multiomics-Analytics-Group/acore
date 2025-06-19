import kmapper as km
import numpy as np
from sklearn import cluster, ensemble


def run_mapper(
    data,
    lenses=["l2norm"],
    n_cubes=15,
    overlap=0.5,
    n_clusters=3,
    linkage="complete",
    affinity="correlation",
):
    """

    :param data:
    :param lenses:
    :param n_cubes:
    :param overlap:
    :param n_clusters:
    :param linkage:
    :param affinity:
    :return:

    """

    X = data._get_numeric_data()
    labels = {i: data.index[i] for i in range(len(data.index))}

    model = ensemble.IsolationForest(random_state=1729)
    model.fit(X)
    lens1 = model.decision_function(X).reshape((X.shape[0], 1))

    # Create another 1-D lens with L2-norm
    mapper = km.KeplerMapper(verbose=0)
    lens2 = mapper.fit_transform(X, projection=lenses[0])

    # Combine both lenses to get a 2-D [Isolation Forest, L^2-Norm] lens
    lens = np.c_[lens1, lens2]

    # Define the simplicial complex
    simplicial_complex = mapper.map(
        lens,
        X,
        clusterer=cluster.AgglomerativeClustering(
            n_clusters=n_clusters, linkage=linkage, affinity=affinity
        ),
    )

    return simplicial_complex, {"labels": labels}
