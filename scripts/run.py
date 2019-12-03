"""
Advanced Network Planning in Four Dimensions - anp4d

Written by Edward Oughton
November 2019
Oxford, UK

"""
import os
import sys
import configparser
import csv

# import math
import fiona
from shapely.geometry import Point, LineString, mapping
import numpy as np
# from random import choice

from rtree import index

from collections import OrderedDict

from anp4d.anp4d import estimate_link_budget

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


def extended_hata(frequency, distance, ant_height, ue_height,
                above_roof, settlement_type, seed_value, iterations):
    """
    Implements the Extended Hata path loss model.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in MHz.
    distance : int
        Distance (d) between transmitter and receiver (kilometres).
    ant_height : int
        Transmitter antenna height (h1) (m, above ground).
    ue_height : int
        Receiver antenna height (h2) (m, above ground).
    above_roof : int
        Whether the path is above or below the roof line (0=below, 1=above).
    settlement_type : string
        General environment (urban/suburban/rural).
    seed_value : int
        Set the seed for the pseudo random number generator
        allowing reproducible stochastic restsults.
    iterations : string
        Specify the number of random numbers to be generated.
        The mean value will be used.

    Returns
    -------
    path_loss : float
        Estimated path loss (dB).

    """
    #find smallest height value
    hm = min(ant_height, ue_height)

    #find largest height value
    hb = max(ant_height, ue_height)

    alpha_hm = (1.1*np.log10(frequency) - 0.7) * min(10, hm) - \
        (1.56*np.log10(frequency) - 0.8) + max(0, (20*np.log10(hm/10)))

    beta_hb = min(0, (20*np.log10(hb/30)))

    if distance <= 20: #units : km
        alpha_exponent = 1

    elif 20 < distance < 100: #units : km
        alpha_exponent = 1 + (0.14 + 1.87e-4 * frequency + \
            1.07e-3 * hb)*(np.log10(distance/20))**0.8
    else:
        raise ValueError('Distance over 100km not compliant')

    ###PART 1####
    #Determine initial path loss based on distance, frequency and environment.
    if distance < 0.04:
        path_loss = (32.4 + (20*np.log10(frequency)) +
            (10*np.log10((distance**2) + ((hb - hm)**2) / (10**6))))

    elif distance >= 0.1:

        if 30 < frequency <= 150:
            path_loss = (69.6 + 26.2*np.log10(150) - 20*np.log10(150/frequency) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                np.log10(distance)**alpha_exponent - alpha_hm - beta_hb)

        elif 150 < frequency <= 1500:
            path_loss = (69.6 + 26.2*np.log10(frequency) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                ((np.log10(distance))**alpha_exponent) - alpha_hm - beta_hb)

        elif 1500 < frequency <= 2000:
            path_loss = (46.3 + 33.9*np.log10(frequency) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                (np.log10(distance))**alpha_exponent - alpha_hm - beta_hb)

        elif 2000 < frequency <= 4000:
            path_loss = (46.3 + 33.9*np.log10(2000) +
                10*np.log10(frequency/2000) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                (np.log10(distance))**alpha_exponent - alpha_hm - beta_hb)

        else:
            raise ValueError('Carrier frequency incorrect for Extended Hata')

        if settlement_type == 'suburban':
            path_loss = (path_loss - 2 * \
                (np.log10((min(max(150, frequency), 2000)/28)))**2 - 5.4)

        elif settlement_type == 'rural': #also called 'open area'
            path_loss = (path_loss - 4.78 * \
                (np.log10(min(max(150, frequency), 2000)))**2 + 18.33 * \
                    np.log10(min(max(150, frequency), 2000)) - 40.94)
        else:
            pass

    elif 0.04 <= distance < 0.1:

        #distance pre-set at 0.1
        l_fixed_distance_upper = (32.4 + (20*np.log10(frequency)) +
              (10*np.log10(0.1**2 + (hb - hm)**2 / 10**6)))

        #distance pre-set at 0.04
        l_fixed_distance_lower = (32.4 + (20*np.log10(frequency)) +
              (10*np.log10(0.04**2 + (hb - hm)**2 / 10**6)))

        path_loss = (l_fixed_distance_lower +
             (np.log10(distance) - np.log10(0.04)) / \
            (np.log10(0.1) - np.log10(0.04)) *
            (l_fixed_distance_upper - l_fixed_distance_lower))

    else:
        # print(distance)
        raise ValueError('Distance over 100km not compliant')

    ###PART 2####
    #determine variation in path loss using stochastic component
    if distance <= 0.04:

        path_loss = path_loss + generate_log_normal_dist_value(frequency,
                                    1, 3.5, iterations, seed_value)

    elif 0.04 < distance <= 0.1:

        if above_roof == 1:
            sigma = (3.5 + ((12-3.5)/0.1-0.04) * (distance - 0.04))
            random_quantity = generate_log_normal_dist_value(frequency,
                                    1, sigma, iterations, seed_value)
            path_loss = (path_loss + random_quantity)

        elif above_roof == 0:
            sigma = (3.5 + ((17-3.5)/0.1-0.04) * (distance - 0.04))
            random_quantity = generate_log_normal_dist_value(frequency,
                                    1, sigma, iterations, seed_value)
            path_loss = (path_loss + random_quantity)

        else:
            raise ValueError('Could not determine if above or below roof line')

    elif 0.1 < distance <= 0.2:

        if above_roof == 1:
            random_quantity = generate_log_normal_dist_value(frequency,
                                1, 12, iterations, seed_value)
            path_loss = (path_loss + random_quantity)
        elif above_roof == 0:
            random_quantity = generate_log_normal_dist_value(frequency,
                                1, 17, iterations, seed_value)
            path_loss = (path_loss + random_quantity)
        else:
            raise ValueError('Could not determine if above or below roof line')

    elif 0.2 < distance <= 0.6:

        if above_roof == 1:
            sigma = (12 + ((9-12)/0.6-0.2) * (distance - 0.02))
            random_quantity = generate_log_normal_dist_value(frequency,
                                1, sigma, iterations, seed_value)
            path_loss = (path_loss + random_quantity)

        elif above_roof == 0:
            sigma = (17 + (9-17) / (0.6-0.2) * (distance - 0.02))
            random_quantity = generate_log_normal_dist_value(frequency,
                                1, sigma, iterations, seed_value)
            path_loss = (path_loss + random_quantity)
        else:
            raise ValueError('Could not determine if above or below roof line')

    elif 0.6 < distance:

        random_quantity = generate_log_normal_dist_value(frequency,
                                1, 12, iterations, seed_value)

        path_loss = (path_loss + random_quantity)

    return round(path_loss, 2)


def generate_log_normal_dist_value(frequency, mu, sigma, iterations, seed_value):
    """
    Generates random values using a lognormal distribution,
    given a specific mean (mu) and standard deviation (sigma).

    https://stackoverflow.com/questions/51609299/python-np-lognormal-gives-infinite-
    results-for-big-average-and-st-dev

    The parameters mu and sigma in np.random.lognormal are not the mean and STD of
    the lognormal distribution. They are the mean and STD of the underlying normal
    distribution.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in MHz.
    mu : int
        Mean of the desired distribution.
    sigma : int
        Standard deviation of the desired distribution.
    iterations : string
        Specify the number of random numbers to be generated.
        The mean value will be used.
    seed_value : int
        Set the seed for the pseudo random number generator
        allowing reproducible stochastic restsults.

    """
    if seed_value == None:
        pass
    else:
        frequency_seed_value = seed_value * frequency * 100

        np.random.seed(int(str(frequency_seed_value)[:2]))

    normal_std = np.sqrt(np.log10(1 + (sigma/mu)**2))
    normal_mean = np.log10(mu) - normal_std**2 / 2

    hs = np.random.lognormal(normal_mean, normal_std, iterations)

    return round(np.mean(hs),2)


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

    modulation_and_coding_lut =[
        # CQI Index, Modulation, Coding rate, Spectral efficiency (bps/Hz), SINR (dB)
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
    path = os.path.join('data','national','fullNetworkWithEdgeIDs.shp')
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
                    capacity_km2 = estimate_link_budget(road_geom, site, frequency, bandwidth,
                        settlement_type, seed_value, iterations, modulation_and_coding_lut)

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
