"""
Path loss calculator
Author: Edward Oughton
Date: November 2019

"""
import numpy as np
from math import pi, sqrt

def path_loss_calculator(model, frequency, distance, ant_height, ant_type,
    building_height, street_width, settlement_type, type_of_sight,
    ue_height, above_roof, indoor, seed_value, iterations):
    """
    Calculate the correct path loss given a range of critera.

    Parameters
    ----------
    model: string
        Specifies which propagation model to use.
    frequency : float
        Carrier band (f) required in MHz.
    distance : float
        Distance between the transmitter and receiver
        in meters.
    ant_height:
        Height of the antenna.
    ant_type : string
        Indicates the type of site antenna (hotspot,
        micro, macro).
    building_height : int
        Height of surrounding buildings in meters (m).
    street_width : float
        Width of street in meters (m).
    settlement_type : string
        Gives the type of settlement (urban, suburban
        or rural).
    type_of_sight : string
        Indicates whether the path is (Non) Line of Sight
        (LOS or NLOS).
    ue_height : float
        Height of the User Equipment.
    above_roof : int
        Indicates if the propagation line is above or below
        building roofs. Above = 1, below = 0.
    indoor : binary
        Indicates if the user is indoor (True) or
        outdoor (False).
    seed_value : int
        Dictates repeatable random number generation.
    iterations : int
        Specifies how many iterations a calculation should
        be run for.

    Returns
    -------
    path_loss : float
        Path loss in decibels (dB)
    model : string
        Type of model used for path loss estimation.

    """
    if model == 'etsi_tr_138_901':
        if 500 < frequency <= 100000:

            path_loss = etsi_tr_138_901(frequency, distance,
                ant_height, ant_type, building_height,
                street_width, settlement_type, type_of_sight,
                ue_height, above_roof, indoor, seed_value,
                iterations
            )

            path_loss = path_loss + outdoor_to_indoor_path_loss(
                    frequency, indoor, seed_value
                )

        else:

            raise ValueError (
                "frequency of {} NOT within correct range".format(
                    frequency
                )
            )

    if model == 'extended_hata':

        path_loss = extended_hata(frequency, distance/1e3, ant_height,
            ue_height, above_roof, settlement_type, seed_value,
            iterations)

    return round(path_loss)


def etsi_tr_138_901(frequency, distance, ant_height, ant_type,
    building_height, street_width, settlement_type, type_of_sight,
    ue_height, above_roof, indoor, seed_value, iterations):
    """
    This is the ETSI 138.901 for the 3GPP TR 38.901 (see the 5G
    study on channel models for frequencies from 0.5 to 100 GHz
    for more information).

    Parameters
    ----------
    frequency : float
        Carrier band (f) required in MHz.
    distance : float
        Distance between the transmitter and receiver (m).
    ant_height:
        Height of the antenna (m).
    ant_type : string
        Indicates the type of site antenna (hotspot,
        micro, macro).
    building_height : int
        Height of surrounding buildings in meters (m).
    street_width : float
        Width of street in meters (m).
    settlement_type : string
        Gives the type of settlement (urban, suburban
        or rural).
    type_of_sight : string
        Indicates whether the path is (Non) Line of Sight
        (LOS or NLOS).
    ue_height : float
        Height of the User Equipment.
    above_roof : int
        Indicates if the propagation line is above or
        below building roofs. Above = 1, below = 0.
    indoor : binary
        Indicates if the user is indoor (True) or
        outdoor (False).
    seed_value : int
        Dictates repeatable random number generation.
    iterations : int
        Specifies how many iterations a calculation
        should be run for.

    Returns
    -------
    path_loss : float
        Path loss in decibels (dB)

    """
    #frequency needs to be in GHz
    fc = frequency / 1e3
    #define speed of light
    c = 3e8

    #calculate breakpoint distances
    #define effective environment height:
    #UMa he = 1 in the documentation.
    he = 1
    #define basestation antenna height
    hbs = ant_height
    #define user equipment antenna height
    hut = ue_height
    #effective basestation antenna height
    h_apost_bs = ant_height - ue_height
    #effective user equipment antenna height
    h_apost_ut = ue_height - he
    #define mean street width
    w = street_width
    #define mean building height
    h = building_height

    #see distance definitions in documentation
    dbp = 2 * pi * hbs * hut * (fc * 1e9) / c
    d_apost_bp = 4 * h_apost_bs * h_apost_ut * (fc*1e9) / c
    d2d_in = 10 #mean d2d_in value
    d2d_out = distance - d2d_in
    d2d = d2d_out + d2d_in
    d3d = sqrt((d2d_out + d2d_in)**2 + (hbs - hut)**2)

    #make sure parameters comply
    check_3gpp_applicability(building_height, street_width,
        ant_height, ue_height
    )

    if settlement_type == 'suburban' or settlement_type == 'rural':

        pl1 = round(
            20*np.log10(40*pi*d3d*fc/3) +
            min(0.03*h**1.72,10) *
            np.log10(d3d) - min(0.044*h**1.72,14.77) +
            0.002*np.log10(h)*d3d +
            generate_log_normal_dist_value(
                fc, 1, 4, iterations, seed_value
            )
        )

        if 10 <= d2d <= dbp:
            if type_of_sight == 'los':
                 return pl1
        pl_rma_los = pl1

        pl2 = round(
            20*np.log10(40*pi*dbp*fc/3) +
            min(0.03*h**1.72,10) *
            np.log10(dbp) - min(0.044*h**1.72,14.77) +
            0.002*np.log10(h)*dbp +
            generate_log_normal_dist_value(
                fc, 1, 4, iterations, seed_value
            ) +
            40*np.log10(d3d / dbp) +
            generate_log_normal_dist_value(
                fc, 1, 6, iterations, seed_value
            )
        )

        if dbp <= d2d <= 10000:
            if type_of_sight == 'los':
                return pl2
        pl_rma_los = pl2

        if type_of_sight == 'nlos':

            pl_apostrophe_rma_nlos = round(
                161.04 - 7.1 * np.log10(w)+7.5*np.log10(h) -
                (24.37 - 3.7 * (h/hbs)**2)*np.log10(hbs) +
                (43.42 - 3.1*np.log10(hbs))*(np.log10(d3d)-3) +
                20*np.log10(fc) -
                (3.2 * (np.log10(11.75*hut))**2 - 4.97) +
                generate_log_normal_dist_value(
                    fc, 1, 8, iterations, seed_value
                )
            )

            pl_rma_nlos = max(pl_apostrophe_rma_nlos, pl_rma_los)

            return pl_rma_nlos

        if d2d > 10000:
            return uma_nlos_optional(frequency, distance,
                ant_height, ue_height, seed_value, iterations)

    elif settlement_type == 'urban':

        pl1 = round(
            28 + 22 * np.log10(d3d) + 20 * np.log10(fc) +
            generate_log_normal_dist_value(
                fc, 1, 4, iterations, seed_value
            )
        )

        if 10 <= d2d <= d_apost_bp:
            if type_of_sight == 'los':
                return pl1
            pl_uma_los = pl1

        pl2 = round(
            28 + 40*np.log10(d3d) + 20 * np.log10(fc) -
            9*np.log10((d_apost_bp)**2 + (hbs-hut)**2) +
            generate_log_normal_dist_value(
                fc, 1, 4, iterations, seed_value
            )
        )

        if d_apost_bp <= d2d <= 5000:
            if type_of_sight == 'los':
                return pl2

        pl_uma_los = pl2

        if type_of_sight == 'nlos':

            if d2d <= 5000:

                pl_apostrophe_uma_nlos = round(
                    13.54 + 39.08 * np.log10(d3d) + 20 *
                    np.log10(fc) - 0.6 * (hut - 1.5) +
                    generate_log_normal_dist_value(
                        fc, 1, 6, iterations, seed_value
                    )
                )

            if d2d > 5000:

                pl_apostrophe_uma_nlos = uma_nlos_optional(
                    frequency, distance, ant_height,
                    ue_height, seed_value, iterations
                )

            pl_uma_nlos = max(pl_apostrophe_uma_nlos, pl_uma_los)

            return pl_uma_nlos

        else:
            return pl_uma_los

    else:
        raise ValueError('Did not recognise settlement_type')

    return print('complete')


def uma_nlos_optional(frequency, distance, ant_height, ue_height,
    seed_value, iterations):
    """
    UMa NLOS / Optional from ETSI TR 138.901 / 3GPP TR 38.901.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in GHz.
    distance : int
        Distance (d) between transmitter and receiver (km).
    ant_height : int
        Transmitter antenna height (h1) (m, above ground).
    ue_height : int
        Receiver antenna height (h2) (m, above ground).
    sigma : int
        Variation in path loss (dB) which is 2.5dB for free space.
    seed_value : int
        Dictates repeatable random number generation.
    iterations : int
        Specifies iterations for a specific calculation.

    Returns
    -------
    path_loss : float
        Path loss in decibels (dB).

    """
    fc = frequency

    d3d = sqrt((distance)**2 + (ant_height - ue_height)**2)

    path_loss = 32.4 + 20*np.log10(fc) + 30*np.log10(d3d)

    random_variation = generate_log_normal_dist_value(
        frequency, 1, 7.8, iterations, seed_value
    )

    return round(path_loss + random_variation)


def check_3gpp_applicability(building_height, street_width,
    ant_height, ue_height):
    """
    Checks that the parameters given conform to the 3gpp model
    assumptions.

    Parameters
    ----------
    building_height : int
        Height of surrounding buildings in meters (m).
    street_width : float
        Width of street in meters (m).
    ant_height:
        Height of the antenna.
    ue_height : float
        Height of the User Equipment.

    Returns
    -------
    overall_compliant : string
        Indicates whether parameters comply (True) or not (False).

    """
    if 5 <= building_height < 50 :
        building_height_compliant = True
    else:
        building_height_compliant = False
        print('building_height not compliant')

    if 5 <= street_width < 50:
        street_width_compliant = True
    else:
        street_width_compliant = False
        print('Street_width not compliant')

    if 10 <= ant_height < 150:
        ant_height_compliant = True
    else:
        ant_height_compliant = False
        print('ant_height not compliant')

    if 1 <= ue_height < 10:
        ue_height_compliant = True
    else:
        ue_height_compliant = False
        print('ue_height not compliant')

    if (building_height_compliant + street_width_compliant +
        ant_height_compliant + ue_height_compliant) == 4:
        overall_compliant = True
    else:
        overall_compliant = False

    return overall_compliant


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
        Whether the path is above or below the roof line
        (0=below, 1=above).
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

    alpha_hm = ((1.1*np.log10(frequency) - 0.7) *
        min(10, hm) - (1.56*np.log10(frequency) - 0.8) +
        max(0, (20*np.log10(hm/10))))

    beta_hb = min(0, (20*np.log10(hb/30)))

    if distance <= 20: #units : km
        alpha_exponent = 1

    elif 20 < distance < 100: #units : km
        alpha_exponent = 1 + (0.14 + 1.87e-4 * frequency + \
            1.07e-3 * hb)*(np.log10(distance/20))**0.8
    else:
        raise ValueError('Distance over 100km not compliant')

    ###PART 1####
    #Determine initial path loss based on distance,
    # frequency and environment.
    if distance < 0.04:
        path_loss = ((32.4 + (20*np.log10(frequency)) +
            (10*np.log10((distance**2) +
            ((hb - hm)**2) / (10**6)))))

    elif distance >= 0.1:

        if 30 < frequency <= 150:
            path_loss = (69.6 + 26.2*np.log10(150) -
                20*np.log10(150/frequency) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                np.log10(distance)**alpha_exponent -
                alpha_hm - beta_hb)

        elif 150 < frequency <= 1500:
            path_loss = (69.6 + 26.2 * np.log10(frequency) -
                13.82 * np.log10(max(30, hb)) +
                (44.9 - 6.55 * np.log10(max(30, hb))) *
                ((np.log10(distance))**alpha_exponent) -
                alpha_hm - beta_hb)

        elif 1500 < frequency <= 2000:
            path_loss = (46.3 + 33.9 * np.log10(frequency) -
                13.82 * np.log10(max(30, hb)) +
                (44.9 - 6.55 * np.log10(max(30, hb))) *
                (np.log10(distance)) ** alpha_exponent -
                alpha_hm - beta_hb)

        elif 2000 < frequency <= 4000:
            path_loss = (46.3 + 33.9*np.log10(2000) +
                10*np.log10(frequency/2000) -
                13.82*np.log10(max(30, hb)) +
                (44.9 - 6.55*np.log10(max(30, hb))) *
                (np.log10(distance))**alpha_exponent -
                alpha_hm - beta_hb)

        else:
            raise ValueError('Frequency incorrect for Extended Hata')

        if settlement_type == 'suburban':
            path_loss = (path_loss - 2 *
                (np.log10((min(max(150, frequency), 2000)/28))) **
                2 - 5.4)

        elif settlement_type == 'rural': #also called 'open area'
            path_loss = (path_loss - 4.78 * \
                (np.log10(min(max(150, frequency), 2000)))**2 +
                18.33 * np.log10(min(max(150, frequency), 2000)) -
                40.94)
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
        raise ValueError('Distance over 100km not compliant')

    ###PART 2####
    #determine variation in path loss using stochastic component
    if distance <= 0.04:

        path_loss = path_loss + generate_log_normal_dist_value(
                    frequency, 1, 3.5, iterations, seed_value)

    elif 0.04 < distance <= 0.1:

        if above_roof == 1:

            sigma = (3.5 + ((12-3.5)/0.1-0.04) * (distance - 0.04))

            random_quantity = generate_log_normal_dist_value(
                            frequency, 1, sigma, iterations, seed_value)

            path_loss = (path_loss + random_quantity)

        elif above_roof == 0:

            sigma = (3.5 + ((17-3.5)/0.1-0.04) * (distance - 0.04))

            random_quantity = generate_log_normal_dist_value(
                            frequency, 1, sigma, iterations, seed_value)

            path_loss = (path_loss + random_quantity)

        else:
            raise ValueError(
                'Could not determine if above or below roof line')

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
            raise ValueError(
                'Could not determine if above or below roof line')

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
            raise ValueError(
                'Could not determine if above or below roof line')

    elif 0.6 < distance:

        random_quantity = generate_log_normal_dist_value(frequency,
                                1, 12, iterations, seed_value)

        path_loss = (path_loss + random_quantity)

    return round(path_loss, 2)


def generate_log_normal_dist_value(frequency, mu, sigma, draws,
    seed_value):
    """
    Generates random values using a lognormal distribution,
    given a specific mean (mu) and standard deviation (sigma).
    https://stackoverflow.com/questions/51609299/
    The parameters mu and sigma in np.random.lognormal are not the mean
    and STD of the lognormal distribution. They are the mean and STD
    of the underlying normal distribution.

    Parameters
    ----------
    frequency : int
        The frequency of the carrier frequency.
    mu : int
        Mean of the desired distribution.
    sigma : int
        Standard deviation of the desired distribution.
    draws : int
        Number of required values.
    seed_value : int
        Dictates repeatable random number generation.

    Returns
    -------
    random_variation : float
        Mean of the random variation over the specified itations.

    """
    if seed_value == None:
        pass
    else:
        frequency_seed_value = seed_value * frequency * 100

        np.random.seed(int(str(frequency_seed_value)[:2]))

    normal_std = np.sqrt(np.log10(1 + (sigma/mu)**2))
    normal_mean = np.log10(mu) - normal_std**2 / 2

    hs = np.random.lognormal(normal_mean, normal_std, draws)

    return round(np.mean(hs),2)


def outdoor_to_indoor_path_loss(frequency, indoor, seed_value):
    """
    ITU-R M.1225 suggests building penetration loss for shadow fading
    can be modelled as a log-normal distribution with a mean and
    standard deviation of 12 dB and 8 dB respectively.

    Parameters
    ----------
    frequency : int
        Carrier band (f) required in MHz.
    indoor : binary
        Indicates if the user is indoor (True) or outdoor (False).
    seed_value : int
        Dictates repeatable random number generation.

    Returns
    -------
    path_loss : float
        Outdoor to indoor path loss in decibels (dB)

    """
    if indoor:

        path_loss = generate_log_normal_dist_value(
            frequency, 12, 8, 1, seed_value
        )

    else:

        path_loss = 0

    return path_loss
