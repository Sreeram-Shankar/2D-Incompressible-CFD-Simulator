import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Circle
import matplotlib as mpl
mpl.rcParams["font.family"] = "Times New Roman"

#function that formats a number to sig figs
def format_sigfigs(x, sigfigs):
    if x == 0: return f"0.{'0'*(sigfigs-1)}"
    else: return f"{x:.{sigfigs}g}"

#function that craetes a 2d scalar plot vs time
def plot_scalar_vs_time(time, quantity, filename, label, ctk_theme):
    #customizes colors based on color mode
    if ctk_theme == "dark":
        bg_color = "#0b1a2e"
        fg_color = "#d8ecff"
        grid_color = "#333333" 
        line_color = "#4faaff"
    else:
        bg_color = "#e6f2ff" 
        fg_color = "#173a5e"   
        grid_color = "#aac4e0" 
        line_color = "#0060cc"

    #creates the plot and checks for nans
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)
    if np.isnan(quantity).all():
        print(f"Skipping graph for {label} — all values are NaN.")
        return
    time = np.array(time)
    quantity = np.array(quantity)
    valid_mask = ~np.isnan(quantity)
    time = time[valid_mask]
    quantity = quantity[valid_mask]
    ax.plot(time, quantity, color=line_color, linewidth=2)

    #configures labels
    ax.set_xlabel("Time", fontsize=12, color=fg_color)
    ax.set_ylabel(label, fontsize=12, color=fg_color)
    ax.tick_params(axis="both", colors=fg_color, labelsize=10)
    ax.set_title(f"{label} vs Time", fontsize=14, color=fg_color)

    #sets the y limits and grid
    ymin = np.min(quantity)
    ymax = np.max(quantity)
    if ymin >= 0:
        y_min = ymin * 0.92
        y_max = ymax * 1.1
    elif ymax <= 0:
        y_min = ymin * 1.1
        y_max = ymax * 0.92
    else:
        y_min = ymin * 1.1
        y_max = ymax * 1.1
    ax.set_ylim(y_min, y_max)
    ax.grid(True, color=grid_color, linestyle="--", linewidth=0.5)
    plt.tight_layout()

    #saves and returns the file
    os.makedirs("graphs", exist_ok=True)
    full_path = os.path.join("graphs", filename)
    fig.savefig(full_path, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return full_path

#function that creates a heatmap
def plot_heatmap_field(field, filename, label, ctk_theme, time, cyl_x=0, cyl_y=0, cyl_r=0.0, L=1, N=2):
    #customizes colors based on color mode
    if ctk_theme == "dark":
        bg_color = "#0b1a2e"
        fg_color = "#d8ecff"
        grid_color = "#333333"
    else:
        bg_color = "#e6f2ff" 
        fg_color = "#173a5e"   
        grid_color = "#aac4e0" 

    #looks for undefined or nan fields
    if np.isnan(field).all() or np.isinf(field).all():
        field = np.full_like(field, np.nan)
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")
    else:
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")

    #plots the heatmap
    fig, ax = plt.subplots(figsize=(10, 4), dpi=150)
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)
    cax = ax.imshow(field, cmap=cmap, origin="lower", interpolation="bicubic")

    #creates a hole for the optional cylinder
    if cyl_r > 0.0:
        radius_pixels = cyl_r / L * (N - 1)
        center_x_pixels = cyl_x / L * (N - 1)
        center_y_pixels = cyl_y / L * (N - 1)
        circle = Circle((center_x_pixels, center_y_pixels), radius_pixels, edgecolor=grid_color, facecolor=grid_color, lw=1.5, zorder=10)
        ax.add_patch(circle)

    #configures labgels and colors
    ax.set_title(f"{label} at time {format_sigfigs(time, 4)}", fontsize=14, color=fg_color)
    ax.tick_params(colors=fg_color)
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label, color=fg_color)
    cbar.ax.yaxis.set_tick_params(color=fg_color)
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=fg_color)
    plt.tight_layout()

    #saves the file and returns the filepath
    os.makedirs("graphs", exist_ok=True)
    full_path = os.path.join("graphs", filename)
    fig.savefig(full_path, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return full_path

#function that creates a streamline plot
def plot_streamline_field(u, v, background_field, filename, label, ctk_theme, time, cyl_x=0, cyl_y=0, cyl_r=0.0, L=1.0, N=2):
    #customizes colors based on color mode
    if ctk_theme == "dark":
        bg_color = "#0b1a2e"
        fg_color = "#d8ecff"
        grid_color = "#333333" 
        stream_color = "#E3E3E3"
    else:
        bg_color = "#e6f2ff" 
        fg_color = "#173a5e"   
        grid_color = "#aac4e0" 
        stream_color = "#292929"

    #looks for undefined or nan fields
    if np.isnan(background_field).all() or np.isinf(background_field).all():
        background_field = np.full_like(background_field, np.nan)
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")
    else:
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")

    #plots the stream plot and the background field
    fig, ax = plt.subplots(figsize=(10, 4), dpi=150)
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)
    cax = ax.imshow(background_field, cmap=cmap, origin="lower", interpolation="bicubic")

    #places the actual streamlines onto the field
    ny, nx = u.shape
    x = np.linspace(0, nx - 1, nx)
    y = np.linspace(0, ny - 1, ny)
    X, Y = np.meshgrid(x, y)
    ax.streamplot(X, Y, u, v, color=stream_color, linewidth=1.0, density=1.2, arrowsize=1.5)

    #creates a hole for the optional cylinder
    if cyl_r > 0.0:
        radius_pixels = cyl_r / L * (N - 1)
        center_x_pixels = cyl_x / L * (N - 1)
        center_y_pixels = cyl_y / L * (N - 1)
        circle = Circle((center_x_pixels, center_y_pixels), radius_pixels, edgecolor=grid_color, facecolor=grid_color, lw=1.5, zorder=10)
        ax.add_patch(circle)

    #configures labgels and colors
    ax.set_title(f"{label} Streamlines at Time {format_sigfigs(time, 4)}", fontsize=14, color=fg_color)
    ax.tick_params(colors=fg_color)
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label, color=fg_color)
    cbar.ax.yaxis.set_tick_params(color=fg_color)
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=fg_color)

    #saves the file and returns the filepath
    plt.tight_layout()
    os.makedirs("graphs", exist_ok=True)
    full_path = os.path.join("graphs", filename)
    fig.savefig(full_path, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return full_path

#function that creates a quiver plot
def plot_quiver_field(u, v, background_field, filename, label, ctk_theme, time, cyl_x=0, cyl_y=0, cyl_r=0.0, L=1.0, N=2, step=4):
    #customizes colors based on color mode
    if ctk_theme == "dark":
        bg_color = "#0b1a2e"
        fg_color = "#d8ecff"
        grid_color = "#333333" 
        arrow_color = "#E3E3E3"
    else:
        bg_color = "#e6f2ff" 
        fg_color = "#173a5e"   
        grid_color = "#aac4e0" 
        arrow_color = "#292929"

    #looks for undefined or nan fields
    if np.isnan(background_field).all() or np.isinf(background_field).all():
        background_field = np.full_like(background_field, np.nan)
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")
    else:
        cmap = plt.cm.get_cmap("viridis").copy()
        cmap.set_bad(color="red")

    #plots the quiver plot and the background field
    fig, ax = plt.subplots(figsize=(10, 4), dpi=150)
    fig.patch.set_facecolor(bg_color)
    ax.set_facecolor(bg_color)
    cax = ax.imshow(background_field, cmap=cmap, origin="lower", interpolation="bicubic")

    #places the actual quivers onto the field
    ny, nx = u.shape
    x = np.arange(0, nx, step)
    y = np.arange(0, ny, step)
    X, Y = np.meshgrid(x, y)
    U = u[::step, ::step]
    V = v[::step, ::step]
    ax.quiver(X, Y, U, V, color=arrow_color, scale=50, headwidth=3)

    #creates a hole for the optional cylinder
    if cyl_r > 0.0:
        radius_pixels = cyl_r / L * (N - 1)
        center_x_pixels = cyl_x / L * (N - 1)
        center_y_pixels = cyl_y / L * (N - 1)
        circle = Circle((center_x_pixels, center_y_pixels), radius_pixels, edgecolor=grid_color, facecolor=grid_color, lw=1.5, zorder=10)
        ax.add_patch(circle)

    #configures labgels and colors
    ax.set_title(f"{label} Quiverplot at Time {format_sigfigs(time, 4)}", fontsize=14, color=fg_color)
    ax.tick_params(colors=fg_color)
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label(label, color=fg_color)
    cbar.ax.yaxis.set_tick_params(color=fg_color)
    plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color=fg_color)

    #saves the file and returns the filepath
    plt.tight_layout()
    os.makedirs("graphs", exist_ok=True)
    full_path = os.path.join("graphs", filename)
    fig.savefig(full_path, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return full_path