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

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']


def plot_map(metric, label, hour, links, flow_max, output_filename):

    ax = links.plot(
        figsize=(10, 12),
        column=metric,
        cmap='YlOrRd',
        norm=matplotlib.colors.Normalize(vmin=0.0001, vmax=flow_max),
        legend=True,
        legend_kwds={'label': label, 'orientation': 'horizontal'}
    )
    plt.title('{:02d}:00_{}'.format(hour, metric), fontsize=16)
    ctx.add_basemap(ax, crs=links.crs)
    plt.savefig(output_filename)
    plt.close()

    return print('Completed {:02d}.png'.format(hour))


def plot_maps(metric, label, flows, flow_max, links, hours):

    for hour, hour_key in enumerate(hours):
        hour_flows = flows[flows.hour == hour_key].copy()
        hour_flows = links.merge(hour_flows, on='road_id_segment')
        plot_name = os.path.join('vis', '{:02d}_{}.png'.format(hour, metric))
        plot_map(metric, label, hour, hour_flows, flow_max, plot_name)

    return print('Plotted all maps')


def generate_gif(metric, gif_path):

    images = []

    filenames = glob.glob(os.path.join(BASE_PATH, '..', 'vis', '*{}.png'.format(metric)))

    for filename in filenames:
        images.append(imageio.imread(filename))

    imageio.mimsave(os.path.join(gif_path), images)

    return print('Generated .gif')


def make_gif(metric, label, path_flow_data, path_link_data, hours, gif_path):

    flows = pd.read_csv(path_flow_data)

    flow_max = flows.capacity_margin.max()

    links = gpd.read_file(path_link_data).rename(columns={'road_id_se': 'road_id_segment'})

    plot_maps(metric, label, flows, flow_max, links, hours)

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

    path_flow_data = os.path.join(BASE_PATH, '..', 'results', 'results.csv')

    path_link_data = os.path.join('data', 'processed', 'chopped_roads.shp')

    plots = [
        ('vehicle_density', 'Vehicle Density (km)'),
        ('capacity_margin', 'Capacity Margin (Mbps/km^2)'),
        ('demand', 'Demand (Mbps/km^2)'),
        ('capacity', 'Capacity (Mbps/km^2)'),
    ]

    for plot in plots:

        metric = plot[0]
        label = plot[1]

        gif_path = os.path.join(BASE_PATH, '..', 'vis', 'movie_{}.gif'.format(metric))

        make_gif(metric, label, path_flow_data, path_link_data, hours, gif_path)

        # pygifsicle.optimize(gif_path)
