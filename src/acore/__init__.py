from importlib.metadata import version

import dsp_pandas  # sets up pandas formatting options

from . import (
    batch_correction,
    decomposition,
    enrichment_analysis,
    exploratory_analysis,
    filter_metabolomics,
    imputation_analysis,
    io,
    multiple_testing,
    normalization,
    types,
    utils,
)

__all__ = [
    "batch_correction",
    "decomposition",
    "dsp_pandas",
    "enrichment_analysis",
    "exploratory_analysis",
    "filter_metabolomics",
    "imputation_analysis",
    "io",
    "multiple_testing",
    "normalization",
    "types",
    "utils",
]
__version__ = version("acore")
