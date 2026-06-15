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
)

__all__ = [
    "dsp_pandas",
    "decomposition",
    "imputation_analysis",
    "io",
    "enrichment_analysis",
    "filter_metabolomics",
    "batch_correction",
    "exploratory_analysis",
    "multiple_testing",
    "normalization",
    "types",
]
__version__ = version("acore")
