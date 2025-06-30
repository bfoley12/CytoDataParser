from .themes import (
    DEFAULT_COLOR_MAPS,
    FACET_STYLE,
    generate_color_map,
    get_color_map
)

# Optionally import your plotting functions too
from .categorical import categorical_plot
from .helpers import rename_legend

__all__ = [
    # Themes
    "DEFAULT_COLOR_MAPS",
    "FACET_STYLE",
    "generate_color_map",
    "get_color_map",

    # Core plotting tools
    "categorical_plot",
    "rename_legend"
]
