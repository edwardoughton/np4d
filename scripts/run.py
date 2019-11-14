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
from itertools import tee
from rtree import index

from collections import OrderedDict

# from anp4d.system_simulator import SimulationManager

CONFIG = configparser.ConfigParser()
CONFIG.read(os.path.join(os.path.dirname(__file__), 'script_config.ini'))
BASE_PATH = CONFIG['file_locations']['base_path']

def get_sites(path):

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

                        # links.append({
                        #     'type': 'Feature',
                        #     'geometry': mapping(line),
                        #     'properties': {
                        #         'edgeid': link
                        #     }
                        # })

                        roads[link_id] = {
                            'road_id': link_id,
                            'geom': line
                        }
                        number += 1
                else:
                    # links.append({
                    #     'type': 'Feature',
                    #     'geometry': mapping(geom),
                    #     'properties': {
                    #         'edgeid': link
                    #     }
                    # })
                    roads[link] = {
                        'road_id': link,
                        'geom': geom
                    }

    return roads


def extended_hata(frequency, distance, ant_height, ue_height, above_roof,
                  settlement_type, seed_value, iterations):
    """
    Implements the Extended Hata path loss model.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in MHz.
    distance : int
        Distance (d) between transmitter and receiver
        required in kilometres.
    ant_height : int
        Transmitter antenna height (h1) (m, above ground).
    ue_height : int
        Receiver antenna height (h2) (m, above ground).
    settlement_type : string
        General environment (urban/suburban/rural)
    seed_value : int
        Set the seed for the pseudo random number generator
        allowing reproducible stochastic restsults.
    iterations : string
        Specify the number of random numbers to be generated.
        The mean value will be used.

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
    #Determine initial path loss according to distance, frequency and environment.
    if distance < 0.04:
        path_loss = (32.4 + (20*np.log10(frequency)) + (10*np.log10((distance**2) +
            ((hb - hm)**2) / (10**6))))

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
            raise ValueError('Could not determine if cell is above or below roof line')

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
            raise ValueError('Could not determine if cell is above or below roof line')

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
            raise ValueError('Could not determine if cell is above or below roof line')

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
    # {'y': 228500.0427, 'x': 432164.8513, 'site_id': 'SP3206428548'},
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


def modulation_scheme_and_coding_rate(sinr, generation, modulation_and_coding_lut):
    """
    Uses the SINR to allocate a modulation scheme and affliated
    coding rate.

    Parameters
    ----------
    sinr : float
        Signal to Interference plus Noise Ratio.
    generation : string
        Generation of cellular technology (e.g. '4G').
    modulation_and_coding_lut : list of tuples
        Lookup table containg sinr and spectral efficiency values.

    """
    for lower, upper in pairwise(modulation_and_coding_lut):
        if lower[0] and upper[0] == generation:

            lower_sinr = lower[5]
            upper_sinr = upper[5]

            if sinr >= lower_sinr and sinr < upper_sinr:
                return lower[4]

            if sinr >= modulation_and_coding_lut[-1][5]:
                return modulation_and_coding_lut[-1][4]

            if sinr < lower_sinr:
                return 0


def pairwise(iterable):
    """
    Return iterable of 2-tuples in a sliding window.

    Parameters
    ----------
    iterable : list
        Sliding window

    Returns
    -------
    list of tuple
        Iterable of 2-tuples
    Example
    -------
    list(pairwise([1,2,3,4]))
        [(1,2),(2,3),(3,4)]

    """
    a, b = tee(iterable)
    next(b, None)

    return zip(a, b)


def estimate_link_budget(road, site, frequency, bandwidth,  settlement_type,
    seed_value, iterations):
    """
    Function for estimating the link budget of a single point.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in MHz.
    bandwidth : int
        Width of the carrier frequency in MHz.
    settlement_type : string
        General environment (urban/suburban/rural)
    seed_value : int
        Set the seed for the pseudo random number generator
        allowing reproducible stochastic restsults.
    iterations : string
        Specify the number of random numbers to be generated.
        The mean value will be used.

    Return
    ------
    mean_capacity_mbps : float
        The average capacity received as Mbps per km^2

    """
    #Get cell site antenna height
    ant_height = 30

    capacity_results = []

    #get user equipment height
    ue_height = 5

    road_geom = road['geom'].interpolate(road['geom'].length / 2)

    #turn path between cell site and user equipment into shapely line object
    line_geom = LineString([(road_geom.x, road_geom.y),(site.x, site.y)])

    # #frequency in MHz, distance in kilometers
    path_loss_dB = extended_hata(frequency, line_geom.length / 1e3,  ant_height, ue_height,
                                0, settlement_type, seed_value, iterations)

    #Equivalent Isotropically Radiated Power (EIRP) - Effective radiated power
    #eirp = site power + site gain - site losses
    eirp = 40 + 16 - 1

    # signal/field strength - received power from the transmitter by a reference antenna
    # at a distance from the transmitting antenna
    #received power = eirp - path_loss - ue_misc_losses + ue_gain - ue_losses
    received_power = eirp - path_loss_dB - 4 + 4 - 4

    #Unwanted in-band interference from other radio antennas
    inteference = -60

    #Unwanted natural, man-made and thermal electromagnetic noise
    #noise parameters
    k = 1.38e-23
    t = 290
    BW = bandwidth*1000000
    noise = 10*np.log10(k*t*1000)+1.5+10*np.log10(BW)

    #calculate the signal to interference plus noise ratio
    sinr = np.log10((10**received_power) / #get the raw linear received power
        ((10**inteference) + #get raw linear sum of interference
            (10**noise))) #get the raw linear noise

    #get the corresponding spectral efficiency achievable with the current sinr
    spectral_efficiency = modulation_scheme_and_coding_rate(
                            sinr, '4G', modulation_and_coding_lut)

    #estimate link budget
    #capacity_mbps = (bits per Hz * channel bandwidth) * 1e6
    link_budget_mbps = (spectral_efficiency * BW) / 1e6

    capacity_results.append(link_budget_mbps)

    mean_capacity_mbps = round(sum(capacity_results) / len(capacity_results))

    return mean_capacity_mbps


def csv_writer(data, directory, filename):
    """
    Write data to a CSV file path
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


def write_shapefile(data, directory, filename, crs):
    """
    Write geojson data to shapefile.

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
    flows, unique_link_ids = load_road_flows(path)[:10]

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

        for flow in flows:
            # print(road_id, flow)
            if int(road_id) == flow['road_id']:
                # {'vehicles': 818, 'road_id': 52722, 'hour': 'FIVEPM'}

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

                    hour = interval_key
                    vehicle_density = flow['vehicles']

                    site = find_closest_site(road, sites)

                    demand_km2 = estimate_demand(vehicle_density, target_capacity, obf)

                    capacity_km2 = estimate_link_budget(road, site, frequency, bandwidth,
                        settlement_type, seed_value, iterations)

                    capacity_margin_km2 = capacity_km2 - demand_km2

                    results.append({
                        'road_id': road_id,
                        'hour': hour,
                        'vehicle_density': vehicle_density,
                        'demand': demand_km2,
                        'capacity': capacity_km2,
                        'capacity_margin': capacity_margin_km2,
                    })

    print('Writing processed sites to .csv')
    csv_writer(results, directory, 'results.csv')

    # print(roads)

    # write_shapefile(roads, directory, 'chopped_roads.shp', crs)
