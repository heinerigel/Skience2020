import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import NullFormatter
import tqdm.notebook as tqdm
import matplotlib.gridspec as _gridspec
from typing import List as _List

def marginal_grid(
    samples,
    dimensions_list: _List[int],
    bins: int = 25,
    show: bool = True,
    colormap_2d=plt.get_cmap("Greys"),
    color_1d="black",
    labels=None,
    limits=None,
):
    number_of_plots = len(dimensions_list)

    plt.figure(figsize=(12, 12))
    gs1 = _gridspec.GridSpec(number_of_plots, number_of_plots)
    gs1.update(wspace=0.05, hspace=0.05)  # set the spacing between axes.

    # Get extent of every set
    dim_range = []
    for i_dim in range(number_of_plots):
        min = samples[dimensions_list[i_dim], :].min()
        max = samples[dimensions_list[i_dim], :].max()
        dim_range.append((min, max))

    for i_plot in range(number_of_plots):
        axis = plt.subplot(gs1[i_plot + (number_of_plots) * i_plot])

        # Modify axes
        if i_plot != number_of_plots - 1:
            axis.set_xticklabels([])
            axis.tick_params(axis="x", which="both", bottom=False, top=False)
            
        axis.set_yticklabels([])
        axis.tick_params(axis="y", which="both", left=False, right=False)
        
        if labels == None:
            label = f"dimension {dimensions_list[number_of_plots-1]}"
        else:
            label = labels[dimensions_list[number_of_plots-1]]
        axis.set_xlabel(label)
        if labels == None:
            label = f"dimension {dimensions_list[number_of_plots-1]}"
        else:
            label = labels[dimensions_list[number_of_plots-1]]
        axis.set_ylabel(label)

        # Plot histogram on diagonal
        axis.hist(
            samples[dimensions_list[i_plot], :],
            bins=bins,
            density=False,
            range=dim_range[i_plot],
            color=color_1d,
        )
        
        if limits is not None:
            axis.set_xlim(limits[i_plot])

        for j_plot in range(i_plot):
            # print(i_plot, j_plot) # grid indices for lower left
            axis = plt.subplot(gs1[j_plot + (number_of_plots) * i_plot], facecolor=(1, 1, 1))

            # Modify axes
            if i_plot != number_of_plots - 1:
                axis.set_xticklabels([])
                axis.tick_params(axis="x", which="both", bottom=False, top=False)
            else:
                if labels == None:
                    label = f"dimension {dimensions_list[j_plot]}"
                else:
                    label = labels[dimensions_list[j_plot]]
                axis.set_xlabel(label)

            if j_plot != 0:
                axis.set_yticklabels([])
                axis.tick_params(axis="y", which="both", left=False, right=False)
            else:
                if labels == None:
                    label = f"dimension {dimensions_list[i_plot]}"
                else:
                    label = labels[dimensions_list[i_plot]]
                axis.set_ylabel(label)

            # Plot 2d marginals
            axis.hist2d(
                samples[dimensions_list[j_plot], :],
                samples[dimensions_list[i_plot], :],
                bins=bins,
                range=[dim_range[j_plot], dim_range[i_plot]],
                cmap=colormap_2d,
            )
            
            if limits is not None:
                axis.set_xlim(limits[j_plot])
                axis.set_ylim(limits[i_plot])

            # print(i_plot, j_plot) # grid indices for lower left
            axis = plt.subplot(gs1[i_plot + (number_of_plots) * j_plot])

            axis.set_xticklabels([])
            axis.tick_params(axis="x", which="both", bottom=False, top=False)
            axis.set_yticklabels([])
            axis.tick_params(axis="y", which="both", left=False, right=False)

            axis.axis("off")

            correlation = np.corrcoef(
                samples[dimensions_list[j_plot], :], samples[dimensions_list[i_plot], :]
            )[1][0]
            axis.text(
                0.5,
                0.5,
                f"{correlation:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=40 * np.abs(correlation),
                transform=axis.transAxes,
            )

    if show:
        plt.show()
