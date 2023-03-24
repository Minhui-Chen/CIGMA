import re
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

def get_xlabel(arg):
    xlabel = arg
    if arg == 'vc':
        xlabel = 'Vacriance proportion of (CT-G, CT-E)'
    elif arg == 'ss':
        xlabel = 'sample size'
    return( xlabel )

def set_xtickname(fig, ax, arg):
    fig.canvas.draw_idle()
    plt.sca( ax )
    locs, labels = plt.xticks()
    news = []
    for label in labels:
        label = label.get_text()
        try:
            if arg == 'vc':
                label = label.split('_')
                label = f'({label[3]}, {label[4]})'
            elif arg == 'ss':
                label = str(int(float(label)))
            news.append( label )
        except Exception as e:
            print(e)
            print(labels)
            news = labels
            break
    plt.xticks( locs, news )
