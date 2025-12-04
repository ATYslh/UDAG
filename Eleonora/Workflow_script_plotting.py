import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

# Workflow script for plotting spatial maps and annual cycles biases

# For spatial maps:


def plotting_structure(
    title,
    models,
    nrows,
    ncols,
    s1,
    s2,
    yn=0.95,
    ts=14,
    ms=10,
    axes_to_set_off=None,
    extent=None,
):
    """
    Plots maps for climatologies or biases. Grid is regular and mapped into a regular projection.

    Parameters:
        title (str): Title for the entire figure.
        models (list of str): Names of the models to be plotted.
        nrows (int): Number of rows in the subplot grid.
        ncols (int): Number of columns in the subplot grid.
        s1, s2 (float): Figure size dimensions.
        yn (float): Vertical position of the title. Default is 0.95.
        ts (int): Font size for the main title. Default is 14.
        ms (int): Font size for the model subplot titles. Default is 10.
        axes_to_set_off (list of int, optional): Indices of axes to hide. Default is None.
        extent (list of float, optional): Geographic extent for the maps. Default is Mediterranean area.

    Returns:
        fig, ax: The matplotlib figure and axes array.
    """

    if axes_to_set_off is None:
        axes_to_set_off = []

    if extent is None:
        extent = [-13, 33, 30, 68]  # Euro-Mediterranean area

    fig, ax = plt.subplots(
        nrows,
        ncols,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(s1, s2),
        sharex="col",
        sharey="row",
    )
    ax = ax.ravel()

    fig.suptitle(title, fontsize=ts, fontweight="bold", y=yn)

    for i, model in enumerate(models):
        ax[i].set_title(model, fontsize=ms, pad=15)
        ax[i].set_extent(extent)

        row = i // ncols
        col = i % ncols
        is_first_in_row = col == 0
        is_last_row = row == nrows - 1

        gl = ax[i].gridlines(
            draw_labels=True, linewidth=0.5, color="gray", alpha=1, linestyle="-"
        )
        gl.xlabel_style = {"size": 8}
        gl.ylabel_style = {"size": 8}

        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = is_last_row
        gl.left_labels = is_first_in_row

        ax[i].coastlines(lw=0.5)
        ax[i].add_feature(cfeature.BORDERS, linestyle="-", linewidth=0.3)

    # Remove axes not used
    total_plots = nrows * ncols
    for i in range(len(models), total_plots):
        ax[i].set_visible(False)

    # Remove specified axes
    for index in axes_to_set_off:
        if 0 <= index < len(ax):
            ax[index].set_visible(False)

    return fig, ax


# Sometimes I used this colorbar:
def personal_colorbar(n_colors):
    # Define the RGB values for the colors
    blue = np.array([0.0, 0.0, 1.0])  # Pure Blue (start)
    white = np.array([1.0, 1.0, 1.0])  # White (middle)
    orange = np.array([1.0, 0.647, 0.0])  # Orange (transition before red)
    red = np.array([1.0, 0.0, 0.0])  # Pure Red (end)

    # Interpolate to generate colors between blue, white, orange, and red
    blue_white = np.linspace(blue, white, n_colors // 2)
    white_orange = np.linspace(white, orange, n_colors // 4)
    orange_red = np.linspace(orange, red, n_colors // 4)

    # Combine all three gradients into one array
    cmap_rgb = np.vstack((blue_white, white_orange, orange_red))

    # Create a ListedColormap
    cmap = mcolors.ListedColormap(cmap_rgb)

    return cmap


# Example of plotting
title = "Biases of x \n GCMs and ICON-RCMs versus ERA5, 1961-1990"
models = [
    "a) CNRM-ESM2-1,$\mathbf{DJF}$",
    "b) ICON-CNRM,$\mathbf{DJF}$",
    "c) CNRM-ESM2-1,$\mathbf{JJA}$",
    "d) ICON-CNRM,$\mathbf{JJA}$",
    "e) EC-Earth-Veg,$\mathbf{DJF}$",
    "f) ICON-ECEARTH,$\mathbf{DJF}$",
    "g) EC-Earth-Veg,$\mathbf{JJA}$",
    "h) ICON-ECEARTH,$\mathbf{JJA}$",
    "i)  MIROC6,$\mathbf{DJF}$",
    "h) ICON-MIROC6,$\mathbf{DJF}$",
    "i) MIROC6,$\mathbf{JJA}$",
    "j) ICON-MIROC6,$\mathbf{JJA}$",
    "k) MPI-ESM1-2-HR,$\mathbf{DJF}$",
    "l) ICON-MPI,$\mathbf{DJF}$",
    "m) MPI-ESM1-2-HR,$\mathbf{JJA}$",
    "n) ICON-MPI,$\mathbf{JJA}$",
]


# cmap = 'BrBG'
cmap = personal_colorbar(n_colors=20)
vmin = -100
vmax = 100
# This only because I want a specific order of the plots:
axes_gcm = [0, 4, 8, 12, 2, 6, 10, 14]
axes_rcm = [1, 5, 9, 13, 3, 7, 11, 15]
index_gcm = [0, 2, 4, 6, 1, 3, 5, 7]
index_rcm = [0, 4, 6, 2, 1, 5, 7, 3]


fig, ax = plotting_structure(title, 0.98, models, 4, 4, 12, 13, axes_to_set_off=None)


for i, idx in zip(axes_gcm, index_gcm):
    cs = ax[i].pcolormesh(
        lon_gcm, lat_gcm, gcm_mean[idx], vmin=vmin, vmax=vmax, cmap=cmap
    )

for j, jdx in zip(axes_rcm, index_rcm):
    cs = ax[j].pcolormesh(
        lon_gcm, lat_gcm, rcm_mean[jdx], vmin=vmin, vmax=vmax, cmap=cmap
    )


ax[1].vlines(
    x=1.16,  # x-position of the vertical line
    ymin=1.2,  # start of line (y-direction)
    ymax=-5,  # end of line (y-direction)
    color="black",
    linewidth=2,
    clip_on=False,
    transform=ax[1].transAxes,  # optional: use Axes coordinates
)

plt.subplots_adjust(wspace=0.35, hspace=-0.2)

positions = [[1.05, 0.26, 0.01, 0.45]]

cb = add_colorbar(
    fig,
    cs,
    position=positions[0],
    label="%",
    labelsize=12,
    orientation="vertical",
    pad=0.05,
    extend="both",
)
# label='\u0394°C'; \u0394$\mathregular{W/m^2}$

plt.savefig("Udag_x.png", bbox_inches="tight")


# ================================
# For Annual cycles:

# Plotting same model in different periods in same subplots  for RCMs Evaluation minus ERA5. Example with net longwave radiation
title = "Annual cycle biases of x \n RCMs Evaluation vs ERA5"
models = [
    "Germany",
    "Alps",
    "BritishIslands",
    "East-Europe",
    "France",
    "IberianPeninsula",
    "Mediterranean",
    "Mid-Europe",
    "Scandinavia",
]


data_to_model = {
    0: "Alps",
    1: "BritishIslands",
    2: "East-Europe",
    3: "France",
    4: "Germany",
    5: "IberianPeninsula",
    6: "Mediterranean",
    7: "Mid-Europe",
    8: "Scandinavia",
}


unit = "$\mathregular{\u0394W/m^2}$"
ylim_min = -80
ylim_max = 80

line_labels = ["1961-1990", "1991-2020"]

# --- Create figure and axes ---
fig, ax = plt.subplots(nrows=3, ncols=4, figsize=(10, 9), sharex="col", sharey="row")
fig.suptitle(title, fontweight="bold", y=1)
fig.tight_layout()

# Turn off unused axes
ax[2, 1].axis("off")
ax[2, 2].axis("off")
ax[2, 3].axis("off")


def plot_subset(data, idx, ax, color="C0", xlabel=False, ylabel=False, label=None):
    x = range(1, 13)
    y = data[idx].values

    ax.plot(x, y, solid_capstyle="butt", color=color)

    ax.grid(True, linewidth=0.5, color="gray", alpha=0.7, linestyle="--")
    ax.axhline(0, color="black", linestyle="--", linewidth=1)

    # Axis limits and ticks
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_xlim(1, 12)
    ax.set_xticks(x)
    ax.tick_params(axis="both", labelsize=9)

    # Optional labels
    if xlabel:
        ax.set_xlabel("Months", fontsize=10)
        ax.xaxis.set_tick_params(labelbottom=True)
    if ylabel:
        ax.set_ylabel(unit, fontsize=10)


# --- NEW helper to plot both datasets with two colors ---
colors = ["C0", "C1"]


def plot_pair(ax, idx, xlabel=False, ylabel=False):
    for color, dataset, label in zip(colors, [data_6190, data_9120], line_labels):
        plot_subset(
            dataset, idx, ax, color=color, xlabel=xlabel, ylabel=ylabel, label=label
        )

    # Use mapping for correct title
    ax.set_title(data_to_model[idx], fontsize=10)


# --- Plot each region
plot_pair(ax[0, 0], 4, ylabel=True)
plot_pair(ax[0, 1], 0)
plot_pair(ax[0, 2], 1)
plot_pair(ax[0, 3], 2)

plot_pair(ax[1, 0], 3, ylabel=True)
plot_pair(ax[1, 1], 5, xlabel=True, ylabel=False)
plot_pair(ax[1, 2], 6, xlabel=True, ylabel=False)
plot_pair(ax[1, 3], 7, xlabel=True, ylabel=False)

plot_pair(ax[2, 0], 8, xlabel=True, ylabel=True)


# --- Shared legend in ax[2,2] ---
from matplotlib.lines import Line2D

# Create custom legend handles
custom_lines = [Line2D([0], [0], color="C0", lw=2), Line2D([0], [0], color="C1", lw=2)]

ax[2, 2].legend(custom_lines, line_labels, loc="center", fontsize=10, frameon=False)
ax[2, 2].axis("off")

plt.savefig("Biases_AC_x.png", bbox_inches="tight")


# ======
# Plotting different models in same subplots, RCMs minus GCMs
# --- Global settings ---
title = "Annual cycle biases of x \n RCMsHist vs GCMsHist, 1991-2020"
models = [
    "Germany",
    "Alps",
    "BritishIslands",
    "East-Europe",
    "France",
    "IberianPeninsula",
    "Mediterranean",
    "Mid-Europe",
    "Scandinavia",
]

data_to_model = {
    0: "Alps",
    1: "BritishIslands",
    2: "East-Europe",
    3: "France",
    4: "Germany",
    5: "IberianPeninsula",
    6: "Mediterranean",
    7: "Mid-Europe",
    8: "Scandinavia",
}


unit = "$\mathregular{W/m^2}$"
ylim_min = -80
ylim_max = 80

# --- Legend labels for lines (4 lines now) ---
# line_labels = ["IFHISTL01_vs_IAEVALL02-T",
#               "IEHISTL01_vs_IAEVALL02-T",
#               "IRHISTL01_vs_IAEVALL02-T",
#               "IMHISTL01_vs_IAEVALL02-T"]

line_labels = [
    "IFHISTL01_vs_CNRM-ESM2-1",
    "IEHISTL01_vs_EC-Earth3-Veg",
    "IRHISTL01_vs_MIROC6",
    "IMHISTL01_vs_MPI-ESM1-2-HR",
]


# line_labels = ["CNRM-ESM2-1_vs_ERA5",
#               "EC-Earth3-Veg_vs_ERA5",
#               "MIROC6_vs_ERA5",
#               "MPI-ESM1-2-HR_vs_ERA5" ]


# --- Create figure and axes ---
fig, ax = plt.subplots(
    nrows=3,
    ncols=4,
    figsize=(10, 9),
    sharex="col",
    sharey="row",
)
fig.suptitle(title, fontweight="bold", y=1)
fig.tight_layout()

# Turn off unused axes
ax[2, 1].axis("off")
ax[2, 2].axis("off")
ax[2, 3].axis("off")


# --- Plotting function ---
def plot_subset(
    data,
    start_idx,
    end_idx,
    ax,
    xlabel=False,
    ylabel=False,
    return_handles=False,
    title=None,
):
    x = range(1, 13)  # Months 1–12
    handles = []
    for i_line, i in enumerate(data[start_idx:end_idx]):
        # Plot line with color C0, C1, ... automatically
        (line,) = ax.plot(x, i.values, solid_capstyle="butt", label=line_labels[i_line])
        handles.append(line)

    # Grid and zero line
    ax.grid(True, linewidth=0.5, color="gray", alpha=0.7, linestyle="--")
    ax.axhline(0, color="black", linestyle="--", linewidth=1)

    # Axis limits and ticks
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_xlim(1, 12)
    ax.set_xticks(x)
    ax.tick_params(axis="both", labelsize=9)

    # Optional labels
    if xlabel:
        ax.set_xlabel("Months", fontsize=10)
        ax.xaxis.set_tick_params(labelbottom=True)
    if ylabel:
        ax.set_ylabel(unit, fontsize=10)

    if return_handles:
        return handles


# --- Plot each region ---
plot_subset(
    data, 16, 20, ax[0, 0], xlabel=False, ylabel=True
)  # <- This subplot has all 4 lines
plot_subset(data, 0, 4, ax[0, 1], xlabel=False, ylabel=False)
plot_subset(data, 4, 8, ax[0, 2], xlabel=False, ylabel=False)
plot_subset(data, 8, 12, ax[0, 3], xlabel=False, ylabel=False)
plot_subset(data, 12, 16, ax[1, 0], xlabel=False, ylabel=True)
plot_subset(data, 20, 24, ax[1, 1], xlabel=True, ylabel=False)
plot_subset(data, 24, 28, ax[1, 2], xlabel=True, ylabel=False)
plot_subset(data, 28, 32, ax[1, 3], xlabel=True, ylabel=False)
plot_subset(data, 32, len(data), ax[2, 0], xlabel=True, ylabel=True)

# --- Shared legend in ax[2,2] ---
# Get handles from ax[0,0], which has all 4 lines
handles, _ = ax[0, 0].get_legend_handles_labels()
ax[2, 2].legend(
    handles=handles, labels=line_labels, loc="center", fontsize=10, frameon=False
)
ax[2, 2].axis("off")


# List of axes in the same order as your current plotting calls
axes_list = [
    ax[0, 0],
    ax[0, 1],
    ax[0, 2],
    ax[0, 3],
    ax[1, 0],
    ax[1, 1],
    ax[1, 2],
    ax[1, 3],
    ax[2, 0],
]

# Loop over axes and assign titles from models list
for ax_i, model_name in zip(axes_list, models):
    ax_i.set_title(model_name, fontsize=10)


plt.savefig("Biases_AC_x.png", bbox_inches="tight")
