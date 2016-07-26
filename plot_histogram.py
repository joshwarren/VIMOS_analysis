## ==================================================================
## 		Plot the histograms of the fields
## ==================================================================
## warrenj 20150112 Routine to plot (nicely) histgrams of the velocity 
## fields.

import numpy as np # for reading files
import matplotlib.pyplot as plt # used for plotting
import os



def plot_histogram(v_binned, galaxy=None, redshift=None, vmin=None, 
    vmax=None, weights=None, title=None, xaxis=None, save=None):

    if vmin is None:
	vmin = np.min(v_binned)
    if vmax is None:
	vmax = np.max(v_binned)

    fig, ax = plt.subplots()

    hist, bins = np.histogram(v_binned, bins=20, range=[vmin,vmax], weights=weights)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)

    if galaxy is not None:
        plt.text(0.02,0.98, "Galaxy: " + galaxy, color='black',
            verticalalignment='top',transform=ax.transAxes)
        if redshift is not None:
            plt.text(0.02,0.93, "Redshift: " + str(round(redshift,3)), 
                color = 'black',verticalalignment='top',
                transform=ax.transAxes)

    if title is not None:
        plt.title(title, y=1.1)

    if xaxis is not None:
        ax.set_xlabel(xaxis)

    if save is not None:
        if not os.path.exists(os.path.dirname(save)):
            os.makedirs(os.path.dirname(save))
        fig.savefig(save, bbox_inches="tight")


    plt.close()
