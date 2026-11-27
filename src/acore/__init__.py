from importlib.metadata import version

import dsp_pandas  # noqa: F401  # sets up pandas formatting options

from . import (
    batch_correction,
    correlation_analysis,
    decomposition,
    differential_regulation,
    drift_correction,
    enrichment_analysis,
    exploratory_analysis,
    filter_metabolomics,
    imputation_analysis,
    io,
    multiple_testing,
    normalization,
)

__all__ = [
    "batch_correction",
    "correlation_analysis",
    "decomposition",
    "drift_correction",
    "dsp_pandas",
    "differential_regulation",
    "enrichment_analysis",
    "exploratory_analysis",
    "filter_metabolomics",
    "imputation_analysis",
    "io",
    "multiple_testing",
    "normalization",
]
__version__ = version("acore")
