from importlib.metadata import version

import dsp_pandas  # sets up pandas formatting options

from . import decomposition, imputation_analysis

__all__ = ["dsp_pandas", "decomposition", "imputation_analysis"]
__version__ = version("acore")
