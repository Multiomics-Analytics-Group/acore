# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
# ---

# %% [markdown]
# # Input Formats
#
# All omics data is tabular. The `rows` are the samples, and the `columns` are the features. The `rows` and `columns` can be annotated with metadata. The data is stored in a `dataframe` object, which is a 2D labeled data structure with columns of potentially different types.
#
# Features:
#
# - `df.index.name` - name of the index (row labels)
# - `df.columns.name` - name of the columns (column labels)
#
# > MultiIndex objects are not supported as most libarires do not support them. metadata is stored in a separate dataframe.
#
#
# The approach is to keep it simple, just use a `dataframe` object, which can be easily serialized to a human readable format, e.g. a `csv` file. The file-format can be anything and the validity of the data is check when data is used, mainly using pandera or a few custom functions.

# %%
