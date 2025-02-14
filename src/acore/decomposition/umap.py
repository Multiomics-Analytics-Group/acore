"""Run UMAP on DataFrame and return result."""
import pandas as pd
import umap


def run_umap(X_scaled, y, random_state=42) -> pd.DataFrame:
    reducer = umap.UMAP(random_state=random_state, n_jobs=1)
    embedding = reducer.fit_transform(X_scaled)
    embedding = pd.DataFrame(
        embedding, index=X_scaled.index, columns=["UMAP 1", "UMAP 2"]
    ).join(y.astype("category"))
    return embedding
