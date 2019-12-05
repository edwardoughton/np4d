"""
Network Planning in Four Dimensions - np4d

Written by Edward Oughton
November 2019
Oxford, UK

"""
import os
import sys
import configparser
import csv

import fiona
from shapely.geometry import Point, LineString, mapping
import numpy as np
from rtree import index

from collections import OrderedDict

from np4d.np4d import estimate_link_budget

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']


def get_sites(path):
    """
    Load cell site locations.

    Parameters
    ----------
    path : string
        path for cell site data.

    Returns
    -------
    sites : list of dicts
        Contains the site_id, and then an x and y coordinate.

    """
    sites = []

    with open(path, 'r') as source:
        reader = csv.DictReader(source)
        for item in reader:
            sites.append({
                'site_id': item['Sitengr'],
                'x': float(item['X']),
                'y': float(item['Y']),
            })

    return sites


def load_road_flows(path):
    """
    Load hourly road flow data.

    Parameters
    ----------
    path : string
        path for road flow data.

    Returns
    -------
    flows : list of dicts
        Contains the road_id, hour and number of vehicles.
    unique_link_ids : list of dicts
        Contains a set of unique road_ids.

    """
    unique_link_ids = set()

    flows = []

    with open(path, 'r') as source:
        reader = csv.DictReader(source)
        for item in reader:
            unique_link_ids.add(int(item['edgeID']))
            flows.append({
                # 'year': intem['year']
                'road_id': int(item['edgeID']),
                'hour': item['hour'],
                'vehicles': int(item['vehicles']),
            })

    return flows, unique_link_ids


def load_roads(path, unique_link_ids):
    """
    Load road shapes.

    Parameters
    ----------
    path : string
        path for road flow data.
    unique_link_ids : list of dicts
        Contains a set of unique road_ids.

    Returns
    -------
    roads : dict of dicts
        Contains the road_id as the key, and then the value is a dict
        containing the road_id and shapely geom.

    """
    roads = {}

    seen_ids = []

    with fiona.open(path) as source:
        for item in source:
            link = int(item['properties']['EdgeID'])
            if link in unique_link_ids:
                seen_ids.append(link)
                geom = LineString(item['geometry']['coordinates'])
                length = geom.length

                if length >= 1000:
                    iterations = round(length / 1000)
                    number = 1
                    for i in range(1, iterations + 1):

                        link_id = str(str(link) + '_' + str(i))

                        if number < iterations:
                            i = i * 1000
                            x = geom.interpolate(i - 1000)
                            y = geom.interpolate(i)
                            line = LineString([x, y])
                        else:
                            i = i * 1000
                            x = geom.interpolate(i - 1000)
                            y = geom.interpolate(length)
                            line = LineString([x, y])

                        roads[link_id] = {
                            'road_id': link_id,
                            'geom': line
                        }
                        number += 1
                else:

                    roads[link] = {
                        'road_id': link,
                        'geom': geom
                    }

    return roads


def find_closest_site(road, sites):
    """
    Finds the closest cell site.

    Parameters
    ----------
    road : dict
        Contains the road_id and shapely geom.
    sites : list of dicts
        Contains the site_id, and then an x and y coordinate.

    Returns
    -------
    nearest_site : object
        The closest cellular site as a shapely object.

    """
    site_geoms = []

    road_geom = road['geom']
    road_centroid = road_geom.interpolate(road_geom.length / 2)

    for site in sites:
        site_geoms.append(Point(site['x'], site['y']))

    idx = index.Index(
        (i, site.bounds, site)
        for i, site in enumerate(site_geoms)
        )

    nearest_site = list(idx.nearest(road_centroid.bounds, 1, objects='raw'))[0]

    return nearest_site


def estimate_demand(vehicle_density, target_capacity, obf):
    """
    Function to estimate the capacity-demand for each section of road.

    Parameters
    ----------
    vehicle_density : float
        The number of vehicles per 1 kilometer stretch of road.
    target_capacity : int
        Target capacity per vehicle in Mbps.
    obf : int
        Overbooking factor.

    """
    demand = vehicle_density * target_capacity / obf

    return round(demand)


def csv_writer(data, directory, filename):
    """
    Write data to a CSV file path.

    Parameters
    ----------
    data : list of dicts
        Data to be written.
    directory : string
        Path to export folder
    filename : string
        Desired filename.

    """
    # Create path
    if not os.path.exists(directory):
        os.makedirs(directory)

    fieldnames = []
    for name, value in data[0].items():
        fieldnames.append(name)

    with open(os.path.join(directory, filename), 'w') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames, lineterminator = '\n')
        writer.writeheader()
        writer.writerows(data)


def convert_shapes_to_geojson(data):
    """
    Converts shapes to geojson.

    Parameters
    ----------
    data : dict of dicts

    Returns
    -------
    links : list of dicts
        Contains geojson format data.

    """
    links = []

    for key, value in data.items():

        links.append({
            'type': 'Feature',
            'geometry': mapping(value['geom']),
            'properties': {
                'road_id_segment': key
            }
        })

    return links


def write_shapefile(data, directory, filename, crs):
    """
    Write geojson data to shapefile.

    Parameters
    ----------
    data : list of dicts
        Data to be written.
    directory : string
        Path to export folder.
    filename : string
        Desired filename.
    crs : string
        Present coordinate reference system (crs).

    """
    prop_schema = []
    for name, value in data[0]['properties'].items():
        fiona_prop_type = next((
            fiona_type for fiona_type, python_type in \
                fiona.FIELD_TYPES_MAP.items() if \
                python_type == type(value)), None
            )

        prop_schema.append((name, fiona_prop_type))

    sink_driver = 'ESRI Shapefile'
    sink_crs = {'init': crs}
    sink_schema = {
        'geometry': data[0]['geometry']['type'],
        'properties': OrderedDict(prop_schema)
    }

    if not os.path.exists(directory):
        os.makedirs(directory)

    with fiona.open(
        os.path.join(directory, filename), 'w',
        driver=sink_driver, crs=sink_crs, schema=sink_schema) as sink:
        for datum in data:
            sink.write(datum)


if __name__ == '__main__':

    ##propagation model can either be:
    ##'etsi_tr_138_901' or
    ##'extended_hata'
    model = 'etsi_tr_138_901'

    modulation_and_coding_lut =[
        # CQI, Modulation, Coding rate, SE (bps/Hz), SINR (dB)
        ('4G', 1, 'QPSK',	0.0762,	0.1523, -6.7),
        ('4G', 2, 'QPSK',	0.1172,	0.2344, -4.7),
        ('4G', 3, 'QPSK',	0.1885,	0.377, -2.3),
        ('4G', 4, 'QPSK',	0.3008,	0.6016, 0.2),
        ('4G', 5, 'QPSK',	0.4385,	0.877, 2.4),
        ('4G', 6, 'QPSK',	0.5879,	1.1758,	4.3),
        ('4G', 7, '16QAM', 0.3691, 1.4766, 5.9),
        ('4G', 8, '16QAM', 0.4785, 1.9141, 8.1),
        ('4G', 9, '16QAM', 0.6016, 2.4063, 10.3),
        ('4G', 10, '64QAM', 0.4551, 2.7305, 11.7),
        ('4G', 11, '64QAM', 0.5537, 3.3223, 14.1),
        ('4G', 12, '64QAM', 0.6504, 3.9023, 16.3),
        ('4G', 13, '64QAM', 0.7539, 4.5234, 18.7),
        ('4G', 14, '64QAM', 0.8525, 5.1152, 21),
        ('4G', 15, '64QAM', 0.9258, 5.5547, 22.7),
    ]

    crs = 'epsg:27700'
    directory = os.path.join(BASE_PATH, 'processed')
    directory_results = os.path.join(BASE_PATH, '..', 'results')

    print('Importing sites data')
    path = os.path.join('data','oxford_cells.csv')
    sites = get_sites(path)

    print('Importing road flow data')
    path = os.path.join('data','link_use_central_oxford.csv')
    flows, unique_link_ids = load_road_flows(path)

    print('Importing road data')
    path = os.path.join('data','shapes','fullNetworkWithEdgeIDs.shp')
    roads = load_roads(path, unique_link_ids)

    frequency = 800
    bandwidth = 10
    settlement_type = 'urban'
    seed_value = 42
    iterations = 1
    target_capacity = 2
    obf = 50

    results = []

    for key, road in roads.items():
        road_id = str(key).split('_')[0]

        for interval_key, interval_name in [
            ('MIDNIGHT', '00'),
            ('ONEAM', '01'),
            ('TWOAM', '02'),
            ('THREEAM', '03'),
            ('FOURAM', '04'),
            ('FIVEAM', '05'),
            ('SIXAM', '06'),
            ('SEVENAM', '07'),
            ('EIGHTAM', '08'),
            ('NINEAM', '09'),
            ('TENAM', '10'),
            ('ELEVENAM', '11'),
            ('NOON', '12'),
            ('ONEPM', '13'),
            ('TWOPM', '14'),
            ('THREEPM', '15'),
            ('FOURPM', '16'),
            ('FIVEPM', '17'),
            ('SIXPM', '18'),
            ('SEVENPM', '19'),
            ('EIGHTPM', '20'),
            ('NINEPM', '21'),
            ('TENPM', '22'),
            ('ELEVENPM', '23')
            ]:

            for flow in flows:
                if int(road_id) == flow['road_id'] and interval_key == flow['hour']:

                    hour = interval_key
                    vehicle_density = flow['vehicles']

                    #find the most likely cell to serve that segment
                    site = find_closest_site(road, sites)

                    #get the demand on this road segment
                    demand_km2 = estimate_demand(vehicle_density, target_capacity, obf)

                    #find centroid of road segment
                    road_geom = road['geom'].interpolate(road['geom'].length / 2)

                    #estimate the capacity of road segment
                    capacity_km2 = estimate_link_budget(model, road_geom, site,
                        frequency, bandwidth, settlement_type, seed_value, iterations,
                        modulation_and_coding_lut)

                    #find the capacity margin of road segment
                    capacity_margin_km2 = capacity_km2 - demand_km2

                    #record results
                    results.append({
                        'road_id': road_id,
                        'road_id_segment': key,
                        'hour': hour,
                        'vehicle_density': vehicle_density,
                        'demand': demand_km2,
                        'capacity': capacity_km2,
                        'capacity_margin': capacity_margin_km2,
                    })

    print('Writing processed sites to .csv')
    csv_writer(results, directory_results, 'results.csv')

    print('Converting roads to geojson')
    roads_geojson = convert_shapes_to_geojson(roads)

    print('Writing roads to .shp')
    write_shapefile(roads_geojson, directory, 'chopped_roads.shp', crs)
