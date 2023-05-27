import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def mycolors(n: int=10, palette: str='muted', desat: float=None) -> list:
    """
    Return colors in hex format from the sns.color_palette('colorblind').

    Parameters
        n:  number of colors. Colorblind palette has 10 colors, if n > 10, cycle the palette.
        palette:    bright, muted, colorblind, pastel, deep, dark
        desat:  Proportion to desaturate each color by

    Returns:
        list of colors in hex format
    """
    colors = sns.color_palette(palette, desat=desat)
    # change the order of purple and grey
    blue, grey =colors[4], colors[7]
    colors[4] = grey
    colors[7] = blue
    colors = (colors*10)[:n]
    colors = [mpl.colors.to_hex(color) for color in colors]
    return colors

def mymarkers() -> list:
    """
    Return a list of markers for plot
    see https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers
    """
    return ['.', 'x', 'P', 's', 'D', 'd', '1', '*']

def snsbox_get_x(x_categories: int, hue_categories: int=1, width:float =0.8) -> np.ndarray:
    """
    Get x coordinates of boxes in Seaborn boxplot.

    Parameters:
        x_categories:   number of categories on x axis
        hue_categories: number of hue categories
        width:  width of all the elements for one level of the major grouping variable. (from seaborn.boxplot with default .8)

    Returns:
        [category1:hue1, category1:hue2, category2:hue1, category2:hue2]
    """
    #print(x_categories, hue_categories, width)
    element_length = width / hue_categories
    if hue_categories == 1:
        return np.arange(x_categories)
    elif hue_categories % 2 == 0:
        shifts = np.arange(hue_categories / 2) * element_length + element_length / 2
        shifts = list(np.flip(shifts*(-1))) + list(shifts)
    else:
        shifts = (np.arange((hue_categories-1) / 2)+1) * element_length
        shifts = list(np.flip(shifts*(-1))) + [0] + list(shifts)
    shifts = [[x] for x in shifts]
    xs = np.repeat(np.atleast_2d(np.arange(x_categories)), hue_categories, axis=0)+np.array(shifts)
    xs = xs.T.flatten()
    return xs

