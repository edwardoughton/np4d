"""
Network Planning in Four Dimensions - np4d

Written by Edward Oughton
November 2019
Oxford, UK

"""
import os
import glob
import configparser

import contextily as ctx
import geopandas as gpd
import imageio
import matplotlib.pyplot as plt
import matplotlib.colors
import pandas as pd
import pygifsicle
import seaborn as sns

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__),'..','scripts','script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']


def plot_map(metric, legend_label, title, hour, roads, flow_min, flow_max, sites,
    output_filename, metric_min, metric_max):

    fig, ax = plt.subplots(figsize=(8, 10))

    plt.rcParams['savefig.pad_inches'] = 0
    plt.autoscale(tight=True)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
    #         hspace = -5, wspace = 0)

    roads.plot(
        # figsize=(8, 10),
        column=metric,
        cmap='RdYlBu',#'hot',#
        norm=matplotlib.colors.Normalize(vmin=-150, vmax=150),
        legend=True,
        legend_kwds={'label': legend_label, 'orientation': 'horizontal'},
        ax=ax
    )

    sites.plot(
        column = 'Cell Site',
        markersize=10,
        legend=True,
        ax=ax
        )

    # plt.legend(, bbox_transform=ax.transAxes)
    plt.title('{:02d}:00 {}'.format(hour, title), fontsize=16)
    ctx.add_basemap(ax, crs=roads.crs)
    plt.savefig(output_filename, pad_inches=0, bbox_inches='tight')
    plt.close()

    return print('Completed {:02d}.png'.format(hour))


def plot_maps(metric, legend_label, title, flows, flow_min, flow_max, roads, sites, hours):

    metric_max = flows[metric].max()
    metric_min = flows[metric].min()

    for hour, hour_key in enumerate(hours):
        hour_flows = flows[flows.hour == hour_key].copy()
        hour_flows = roads.merge(hour_flows, on='road_id_segment')
        plot_name = os.path.join('vis', 'images', '{:02d}_{}.png'.format(hour, metric))
        plot_map(metric, legend_label, title, hour, hour_flows, flow_min, flow_max, sites,
            plot_name, metric_min, metric_max)

    return print('Plotted all maps')


def generate_gif(metric, gif_path):

    images = []

    filenames = glob.glob(
        os.path.join(BASE_PATH, '..', 'vis', 'images','*{}.png'.format(metric))
    )

    for filename in filenames:
        images.append(imageio.imread(filename))

    imageio.mimsave(os.path.join(gif_path), images)

    return print('Generated .gif')


def make_gif(metric, legend_label, title, path_flows, path_roads, path_sites, hours, gif_path):

    flows = pd.read_csv(path_flows)

    flow_min = flows[metric].min()
    flow_max = flows[metric].max()

    roads = gpd.read_file(path_roads).rename(columns={'road_id_se': 'road_id_segment'})

    sites = gpd.read_file(path_sites)

    plot_maps(metric, legend_label, title, flows, flow_min, flow_max, roads, sites, hours)

    generate_gif(metric, gif_path)

    return print('Complete')


if __name__ == '__main__':

    hours = [
        'MIDNIGHT',
        'ONEAM',
        'TWOAM',
        'THREEAM',
        'FOURAM',
        'FIVEAM',
        'SIXAM',
        'SEVENAM',
        'EIGHTAM',
        'NINEAM',
        'TENAM',
        'ELEVENAM',
        'NOON',
        'ONEPM',
        'TWOPM',
        'THREEPM',
        'FOURPM',
        'FIVEPM',
        'SIXPM',
        'SEVENPM',
        'EIGHTPM',
        'NINEPM',
        'TENPM',
        'ELEVENPM',
    ]

    path_flows = os.path.join(BASE_PATH, '..', 'results', 'results.csv')

    path_roads = os.path.join('data', 'processed', 'chopped_roads.shp')

    path_sites = os.path.join('data', 'processed', 'sites.shp')

    plots = [
        # ('vehicle_density', 'Vehicle Density (km)'),
        ('capacity_margin', 'Capacity Margin (Mbps/km^2)', 'Capacity Margin'),
        # ('demand', 'Demand (Mbps/km^2)'),
        # ('capacity', 'Capacity (Mbps/km^2)'),
    ]

    for plot in plots:

        metric = plot[0]
        legend_label = plot[1]
        title = plot[2]

        gif_path = os.path.join(
            BASE_PATH, '..', 'vis', 'movies', 'movie_{}.gif'.format(metric)
        )

        make_gif(metric, legend_label, title, path_flows, path_roads, path_sites, hours, gif_path)

        pygifsicle.optimize(gif_path)
