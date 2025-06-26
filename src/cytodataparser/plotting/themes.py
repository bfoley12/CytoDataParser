from plotly.colors import qualitative
import itertools

# === STATIC DEFAULTS === #
DEFAULT_COLOR_MAPS = {
    "strain_source": {
        "B6: Recipient": "blue",
        "F1: Recipient": "red",
        "B6: Donor": "darkblue",
        "F1: Donor": "pink"
    },
    "strain": {
        "B6": "blue",
        "F1": "red"
    },
    "source":{
        "Recipient": "blue",
        "Donor": "red"
    },
    "lobe": {
        "L": "blue",
        "Left": "blue",
        "R": "red",
        "Right": "red"
    }
}

FACET_STYLE = {
    "spacing": 0.05,
    "margin": dict(l=40, r=10, t=40, b=40),
    "showticklabels": True
}


def generate_color_map(values, palette=qualitative.Plotly):
    """
    Generates a color map from a list of values using a given Plotly color palette.

    Parameters:
    - values (list): Unique categories to map
    - palette (list): List of colors (default is plotly qualitative)

    Returns:
    - dict: value -> color
    """
    unique_vals = list(dict.fromkeys(values))  # preserve order, dedupe
    color_cycle = itertools.cycle(palette)
    return {val: next(color_cycle) for val in unique_vals}

def get_color_map(field, values=None):
    """
    Gets the appropriate color map. If known field, returns default. Else, generates one.

    Parameters:
    - field (str): Column name
    - values (list): If not a known field, pass values to auto-generate map

    Returns:
    - dict: value -> color
    """
    lower_field = field.lower()
    if lower_field in DEFAULT_COLOR_MAPS:
        return DEFAULT_COLOR_MAPS[lower_field]
    elif values is not None:
        return generate_color_map(values)
    else:
        raise ValueError(f"No default color map for field '{field}', and no values provided.")

__all__ = [
    "DEFAULT_COLOR_MAPS",
    "FACET_STYLE",
    "generate_color_map",
    "get_color_map"
]
