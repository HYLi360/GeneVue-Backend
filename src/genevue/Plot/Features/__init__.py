"""
A ported version of DNA Features Viewer.

This is a brief introduction of code organization from original repository (fixed some typos).

----------
This document walks you through the DNA Features Viewer code. Please request changes if anything is unclear.

- **GraphicFeature.py** implements a class for defining a *GraphicFeature*, which is an annotation (start, end, strand,
  label) with graphical properties (color, line width, font family...).

- **GraphicRecord/** implements the *GraphicRecord* class, which can plot a set of *GraphicFeatures* using Matplotlib or
  Bokeh. To keep file sizes acceptable, many methods are implemented in separate files (*bokeh_plots.py*,
  *matplotlib_plots.py*) and added to *GraphicRecord* via class mixins.

- **CircularGraphicRecord/** implements the *GraphicRecord* class, which inherits from *GraphicRecord* but draws features
  circularly using custom Matplotlib patches called "arrow-wedge" (defined in file *ArrowWedge.py*).

- **compute_features_levels.py** implements the algorithm for deciding the levels on which the different features (and
  annotations) are drawn.

- **biotools.py** implements generic biology-related methods (reverse complement, annotation of Biopython records, etc.).

----------
This version retains the vast majority of the code in the codebase, but the functionality of *biotools.py* has been
replaced by another module, *genevue.Sequences*.
"""

__version__ = "3.1.5"

from .GraphicRecord import GraphicRecord
from .CircularGraphicRecord import CircularGraphicRecord
from .GraphicFeature import GraphicFeature
from .BiopythonTranslator import (
    BiopythonTranslator,
    BlackBoxlessLabelTranslator,
)
from .biotools import load_record, annotate_biopython_record

__all__ = [
    "GraphicRecord",
    "CircularGraphicRecord",
    "GraphicFeature",
    "BiopythonTranslator",
    "BlackBoxlessLabelTranslator",
    "annotate_biopython_record",
    "__version__",
]
