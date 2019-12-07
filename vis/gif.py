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


def plot_map(metric, label, hour, roads, flow_max, sites, output_filename,
    metric_min, metric_max):

    # fig, (a0, a1) = plt.subplots(
    #     2,
    #     1,
    #     gridspec_kw={'height_ratios': [3, 12]},
    #     figsize = (12, 12)
    # )

    # df = pd.DataFrame(roads)

    # # Density Plot and Histogram of all arrival delays
    # a0 = sns.distplot(df[metric],
    #     hist=True,
    #     # kde=True,
    #     bins=10, #int(180/5),
    #     color = 'darkblue',
    #     hist_kws={'edgecolor':'black'},
    #     kde_kws={'linewidth': 4},
    #     ax=a0
    #     )

    # a0.set_ylim(0, 0.08)
    # a0.set_xlim(-100, 100)

    # a1 = roads.plot(
    #     figsize=(10, 12),
    #     column=metric,
    #     cmap='YlOrRd',
    #     norm=matplotlib.colors.Normalize(vmin=0.0001, vmax=flow_max),
    #     legend=True,
    #     legend_kwds={'label': label, 'orientation': 'horizontal'},
    #     ax=a1
    # )
    # a1 = sites.plot(
    #     column = 'Cell Site',
    #     markersize=10,
    #     legend=True,
    #     ax=a1
    #     )

    # plt.subplots_adjust(bottom=0.1)
    # plt.title('{:02d}:00_{}'.format(hour, metric), fontsize=16)
    # ctx.add_basemap(a1, crs=roads.crs)
    # fig.tight_layout()
    # plt.savefig(output_filename)
    # plt.close()

    fig, ax = plt.subplots(figsize = (12, 12))

    roads.plot(
        figsize=(10, 12),
        column=metric,
        cmap='YlOrRd',
        norm=matplotlib.colors.Normalize(vmin=0.0001, vmax=flow_max),
        legend=True,
        legend_kwds={
            'label': label,
            'orientation': 'horizontal'},
        ax=ax
    ).legend(bbox_to_anchor=(-5, -10))

    # plt.subplots_adjust(bottom=0.99)
    # plt.title('{:02d}:00_{}'.format(hour, metric), fontsize=16)
    # ctx.add_basemap(a1, crs=roads.crs)
    # fig.tight_layout()
    # plt.label(bbox_to_anchor = [0.5, 0.2])
    ax.legend(bbox_to_anchor=(4, 1), loc=1, borderaxespad=-2)
    plt.savefig(output_filename, bbox_to_anchor = [-1, 0.2]) #bbox_inches='tight', pad_inches=-0.3
    plt.close()

    return print('Completed {:02d}.png'.format(hour))


def plot_maps(metric, label, flows, flow_max, roads, sites, hours):

    metric_max = flows[metric].max()
    metric_min = flows[metric].min()

    for hour, hour_key in enumerate(hours):
        hour_flows = flows[flows.hour == hour_key].copy()
        hour_flows = roads.merge(hour_flows, on='road_id_segment')
        plot_name = os.path.join('vis', 'images', '{:02d}_{}.png'.format(hour, metric))
        plot_map(metric, label, hour, hour_flows, flow_max, sites, plot_name, metric_min, metric_max)

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


def make_gif(metric, label, path_flows, path_roads, path_sites, hours, gif_path):

    flows = pd.read_csv(path_flows)

    flow_max = flows[metric].max()

    roads = gpd.read_file(path_roads).rename(columns={'road_id_se': 'road_id_segment'})

    sites = gpd.read_file(path_sites)

    plot_maps(metric, label, flows, flow_max, roads, sites, hours)

    generate_gif(metric, gif_path)

    return print('Complete')


if __name__ == '__main__':

    hours = [
        'MIDNIGHT',
        # 'ONEAM',
        # 'TWOAM',
        # 'THREEAM',
        # 'FOURAM',
        # 'FIVEAM',
        # 'SIXAM',
        # 'SEVENAM',
        # 'EIGHTAM',
        # 'NINEAM',
        # 'TENAM',
        # 'ELEVENAM',
        # 'NOON',
        # 'ONEPM',
        # 'TWOPM',
        # 'THREEPM',
        # 'FOURPM',
        # 'FIVEPM',
        # 'SIXPM',
        # 'SEVENPM',
        # 'EIGHTPM',
        # 'NINEPM',
        # 'TENPM',
        # 'ELEVENPM',
    ]

    path_flows = os.path.join(BASE_PATH, '..', 'results', 'results.csv')

    path_roads = os.path.join('data', 'processed', 'chopped_roads.shp')

    path_sites = os.path.join('data', 'processed', 'sites.shp')

    plots = [
        # ('vehicle_density', 'Vehicle Density (km)'),
        ('capacity_margin', 'Capacity Margin (Mbps/km^2)'),
        # ('demand', 'Demand (Mbps/km^2)'),
        # ('capacity', 'Capacity (Mbps/km^2)'),
    ]

    for plot in plots:

        metric = plot[0]
        label = plot[1]

        gif_path = os.path.join(
            BASE_PATH, '..', 'vis', 'movies', 'movie_{}.gif'.format(metric)
        )

        make_gif(metric, label, path_flows, path_roads, path_sites, hours, gif_path)

        # pygifsicle.optimize(gif_path)
