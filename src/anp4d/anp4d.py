"""
Network Planning in Four Dimensions - np4d

Written by Edward Oughton
November 2019
Oxford, UK

"""
from shapely.geometry import Point, LineString, mapping
import numpy as np
from itertools import tee

from anp4d.path_loss import path_loss_calculator


def estimate_link_budget(receiver, site, frequency, bandwidth, settlement_type,
    seed_value, iterations, modulation_and_coding_lut):
    """
    Function for estimating the link budget of a single point.

    Parameters
    ----------
    receiver : Shapely Point object
        The point geometry for the receiver location.
    site : Shapely Point object
        The point geometry for the cell site location.
    frequency : int
        Carrier band (f) required in MHz.
    bandwidth : int
        Width of the carrier frequency in MHz.
    settlement_type : string
        General environment (urban/suburban/rural).
    seed_value : int
        Set the seed for the pseudo random number generator
        allowing reproducible stochastic restsults.
    iterations : string
        Specify the number of random numbers to be generated.
        The mean value will be used.
    modulation_and_coding_lut : list of tuples
        Lookup table containg sinr and spectral efficiency values.

    Return
    ------
    mean_capacity_mbps : float
        The average capacity received as Mbps per km^2

    """
    capacity_results = []

    #turn path between cell site and user equipment into shapely line object
    line_geom = LineString([(receiver.x, receiver.y),(site.x, site.y)])

    frequency = frequency
    distance = line_geom.length
    ant_height = 30
    ant_type = 'macro'
    building_height = 20
    street_width = 20
    settlement_type = 'urban'
    type_of_sight = 'los'
    ue_height = 5
    above_roof = 0
    indoor = 0
    seed_value = 42
    iterations = 20

    # #frequency in MHz, distance in kilometers
    path_loss_dB = path_loss_calculator(frequency, distance, ant_height, ant_type,
        building_height, street_width, settlement_type, type_of_sight,
        ue_height, above_roof, indoor, seed_value, iterations)

    #Equivalent Isotropically Radiated Power (EIRP) - Effective radiated power
    #eirp = site power + site gain - site losses
    eirp = 40 + 16 - 1

    # signal/field strength - received power from the transmitter by a reference antenna
    # at a distance from the transmitting antenna
    #received power = eirp - path_loss - ue_misc_losses + ue_gain - ue_losses
    received_power = eirp - path_loss_dB[0] - 4 + 4 - 4

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
