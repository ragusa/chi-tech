#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 10:34:55 2023

@author: jean.ragusa
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import csv as csv
import matplotlib.pyplot as plt
plt.close('all')

# %%


def read_csv_file(file_path):
    data = []

    try:
        with open(file_path, 'r', newline='') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                if len(row) == 3:
                    try:
                        # Convert the values to floats
                        converted_row = [float(value) for value in row]
                        data.append(converted_row)
                    except ValueError as ve:
                        print(f"Skipping row due to conversion error: {ve}")
                else:
                    print(
                        f"Skipping row with incorrect number of columns: {row}")

        return data

    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


# %%
def plot_angular_quadrature(file_path):
    # Read csv data
    csv_data = read_csv_file(file_path)
    qdata = np.asarray(csv_data)
    # determine if 2d or 3d quadrature data by looking at sum of mu's
    if np.sum(qdata[:, 1]) > 0.01:
        dim = 2
    # get number of directions
    n_dir = qdata.shape[0]
    # compute omega
    omega = np.zeros((n_dir, 3))
    # Omega_z
    mu = np.copy(qdata[:, 1])
    omega[:, 2] = mu[:]
    # Omega_x
    omega[:, 0] = np.sqrt(1-mu**2) * np.cos(qdata[:, 2])
    # Omega_y
    omega[:, 1] = np.sqrt(1-mu**2) * np.sin(qdata[:, 2])

    # create figure with options
    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(111, projection='3d')

    # transparent sphere data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    # Plot the surface
    ax.plot_surface(x, y, z, color='b', alpha=0.1)

    # create axes
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def do_3d_projection(self, renderer=None):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        def draw(self, renderer=None):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)

    # a = Arrow3D([-1, 1], [0, 0], [0, 0], mutation_scale=10,
    #             lw=1, arrowstyle='->', color="chartreuse")
    # ax.add_artist(a)
    # a = Arrow3D([0, 0], [-1, 1], [0, 0], mutation_scale=10,
    #             lw=1, arrowstyle='->', color="indigo")
    # ax.add_artist(a)
    # a = Arrow3D([0, 0], [0, 0], [-1, 1], mutation_scale=10,
    #             lw=1, arrowstyle='->', color="darkorange")
    # ax.add_artist(a)

    # http://zyxue.github.io/2017/04/13/matplotlib-color-marker-combinations.html
    import itertools
    colors = itertools.cycle(["r", "g", "b", "k", "m", "c", "y", "crimson"])
    # https://matplotlib.org/3.1.0/gallery/color/named_colors.html

    # plot quadrature itself
    d_skip = int(0)
    n_octants = int(2**dim)
    n_dir_per_oct = int(n_dir / n_octants)
    print("n octants=", n_octants, "ndir per octant=", n_dir_per_oct)
    for oc in range(n_octants):
        clr = next(colors)
        d_skip = oc*n_dir_per_oct
        for d in range(n_dir_per_oct):
            om = omega[d+d_skip, :]
            ax.plot3D([0, om[0]], [0, om[1]], [
                      0, om[2]], c=clr, linewidth=0.75)

    polar_level = np.unique(mu)
    for p in range(len(polar_level)):
        r = polar_level[p]
        x = np.sqrt(1-r**2)*np.cos(u)
        y = np.sqrt(1-r**2)*np.sin(u)
        z = r*np.ones(len(u))
        ax.plot3D(x, y, z, 'grey', linestyle="dashed", linewidth=0.5)

    ax.view_init(30, 70)
    plt.show()
    print(omega)


# %%
# Usage example
file_path = "qdata_na2_np4.csv"  # Replace with the path to your CSV file

plot_angular_quadrature(file_path)
