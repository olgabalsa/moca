import argparse
import numpy as np
import sys
import os
import scripts.st_to_mass as st_to_mass # type: ignore
import scripts.rv_functions as rv_functions # type: ignore
import scripts.forecaster.mr_forecast as mr # type: ignore
from astropy.time import Time
from termcolor import cprint, colored
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.stats import truncnorm
import pandas as pd


def print_param_planet(i_pl, num_pl, param, simul, name_param, unit_param, mp_from_rp, erp = None, erpsim = None):
    try:
        if not erp:
            erp = np.array([None] * num_pl)
    except:
        pass
    
    # TO DO: revisar estas lineas
    if name_param == 'Mp' and mp_from_rp[i_pl] == True:
        text1 = colored(f'Mp = {np.round(param[i_pl], 2)} {unit_param} -', 'white')
        text2 = colored('computed from Rp', 'cyan')
        print(text1, text2)
    elif erp[i_pl]:
        if erpsim[i_pl]:
            if simul[i_pl]:
                text1 = colored(f'Rp = {np.round(param[i_pl], 2)} +- {np.round(erp[i_pl], 2)} {unit_param} -', 'white')
                text2 = colored('randomly chosen', 'yellow')
                print(text1, text2)
            else:
                text1 = colored(f'Rp = {np.round(param[i_pl], 2)} +- {np.round(erp[i_pl], 2)} {unit_param} -', 'white')
                text2 = colored('used eRp = 0.1 x Rp as default.', 'yellow')
                text3 = colored('NOTE:', 'blue')
                text4 = colored('eRp choice strongly affect the mass estimation', 'white')
                print(text1, text2)
        else:
            print(f'Rp = {np.round(param[i_pl], 2)} +- {np.round(erp[i_pl], 2)} {unit_param}')
    else:
        if simul[i_pl] == 1:
            text1 = colored(f'{name_param}{i_pl+1} = {np.round(param[i_pl], 2)} {unit_param} -', 'white')
            text2 = colored('randomly chosen', 'yellow')
            print(text1, text2)
        else:
            print(f'{name_param}{i_pl+1} = {np.round(param[i_pl], 2)} {unit_param}')


def bool_str(val, param):
    if val == 'False' or val == False:
        return False
    elif val == 'True':
        return True
    elif not val:
        return None
    else:
        cprint(f'ERROR: {param} input must be True or False', 'red')
    

def get():

    #--------------------------------------------
    parser = argparse.ArgumentParser()

    #----------
    # Star info
    #----------
    parser.add_argument("-star", "--star_name", type = str, default = 'Unknown', help = "Name of the star to be resolved with SIMBAD [e.g., -star TOI-969]. Needed to consider the visibility of the target, alternatively to its coordinates (-ra and -dec).")
    parser.add_argument("-ra", "--right_ascension", type = float, default = None, help = "RA of the star in deg units in demical format [e.g., -ra 7.6758]. Needed to consider the visibility of the target, alternatively to its SIMBAD name (-s).")
    parser.add_argument("-dec", "--declination", type = float, default = None, help = "DEC of the star in deg units in demical format [e.g., -dec 2.0989]. Needed to consider the visibility of the target, alternatively to its SIMBAD name (-s).")
    
    parser.add_argument("-ms", "--stellar_mass", type = float, default = None, help = "Stellar mass in solar units [e.g., -ms 0.734].")
    parser.add_argument("-st", "--spectral_type", type = str, default = None, help = "Stellar spectral type for main sequence stars [e.g., -st K4]. Used to get the stellar mass if not provided.")

    parser.add_argument("-act", "--activity_level", type = int, default = None, help = "Level of stellar activity [e.g., -act 1]. Options: from 0 to 3 in increasing activity level (Inactive, Quiet, Moderate and Loud).")
    parser.add_argument("-pr", "--P_rotation", type = float, default = None, help = "Stellar rotational period in days [e.g., -pr 24.6]. This value is used for including the red noise due to stellar activity if -act different from 0.")
    parser.add_argument("-pm", "--P_magnetic_cycle", type = float, default = None, help = "Period of the magnetic cycle in days [e.g., -pm 100]. This value is used for including the red noise due to stellar activity if -act different from 0.")


    #------------
    # Planet info
    #------------
    parser.add_argument("-np", "--number_planets", type = int, default = None, help = "Number of planets [e.g., -np 2].")
    parser.add_argument("-p", "--orbital_period", nargs = '+', default = None, help = "Orbital period of the planet(s) in days [e.g., -p 1.824 1700].")
    parser.add_argument("-t0", "--transiting_time", nargs = '+', default = None, help = "Transiting time of the planet(s) in BJD [e.g., -t0 2459248.37709 nan]. Note that 'nan' can be used in any planetary parameter but the number of planets to use a random value. If none of the -t0 were provided, both would be random.")
    parser.add_argument("-e", "--eccentricity", nargs = '+', default = None, help = "Orbital eccentricity(/ies) [e.g., -e nan 0.628].")
    parser.add_argument("-w", "--arg_periastron", nargs = '+', default = None, help = "Argument(s) of periastron in degrees [e.g., -w nan 208.5].") 
    parser.add_argument("-i", "--inclination", nargs = '+', default = None, help = "Orbital inclination(s) in degrees [e.g., 86.75 nan].")
    
    parser.add_argument("-mp", "--planetary_mass", nargs = '+', default = None, help = "Planetary mass(es) in the units specified in the -u input [e.g., -mp 9.1 3590].")
    parser.add_argument("-rp", "--planetary_radius", nargs = '+', default = None, help = "Planetary radius in the units specified in the -u input [e.g., -rp 2.765 nan]. Used to infer the mass if -mp not given.")
    parser.add_argument("-erp", "--unc_planetary_radius", nargs = '+', default = None, help = "Uncertainty associated to the planetary radius (in the units specified in the -u input) [e.g., -erp 0.097 nan]. Needed to infer the mass. It is set to 0.2 x Rp if not provided.")
    parser.add_argument("-u", "--units", type = str, default = None, help = "Units for the planetary mass and radius [e.g., -u e]. Two options: [e, j] for Earth and Jupiter units, respectively.")


    #-------------------
    # Observational info
    #-------------------
    parser.add_argument("-erv", "--rv_uncertainty", type = float, default = None, help = "Median of the RV uncertainties to be simulated in m/s [e.g., -erv 3.5].")
    parser.add_argument("-obs", "--observatory", type = str, default = None, help = "Name of the observatory [e.g., -obs 'La Silla']. Used to compute the target visibility (if target name or coordinates are also given) and customize the weather loss (if not turned off with -wl False). Options: ['CAHA', 'Large Binocular', 'La Silla', 'Lick', 'McDonald', 'Mauna Kea', 'OHP', 'Paranal', 'Roque Muchachos', 'Whipple']")
    parser.add_argument("-wl", "--weather_loss", type = str, default = None, help = "Boolean to consider weather loss [e.g., -wl False]. If the observatory is also provided, the weather loss is based on its statistics per month, otherwise, a homogeneous loss of 30 percent is considered.")
    parser.add_argument("-cwlp", "--custom_weather_loss_percentage", type = float, default = None, help = "Customized weather loss percentage [e.g., -cwlp 10]. An homogenous weather loss of the given percentage will be considered regardless of the observatory.")
    parser.add_argument("-id", "--init_date", type = str, default = None, help = "Date to start the observations in 'YYYY-MM-DD' format [used to compute the target visibility if observatory and target name/coordinates are provided].")
    parser.add_argument("-ed", "--end_date", type = str, default = None, help = "Date to end the observations in 'YYYY-MM-DD' format [used to compute the target visibility if observatory and target name/coordinates are provided].")
    parser.add_argument("-ma", "--minimum_altitude", type = float, default = None, help = "Minumum altitude of the target (in degrees) to be observed [e.g., -ma 30].")


    #-----------------------
    # Observational strategy
    #-----------------------
    parser.add_argument("-sn", "--strategy_name", type = str, default = None, help = "Strategy name [e.g., -ns mon]. Options: ['mon', 'blo', 'jds'] for monotonic, blocks, or for JDs already scheduled, respectively. For more information on each strategy see the documentation at (TO BE ADDED)") # TD

    # FOR Monotonic & Blocks
    parser.add_argument("-nm", "--number_max_rvs", type = int, default = None, help = "Maximum total number of RV measurements to be simulated [e.g., -nm 30]. Note that the simulated dates can be below this number due to wheather loss and visibility conditions (if those are considered).")
    parser.add_argument("-nn", "--number_per_night_rvs", type = int, default = None, help = "Number of RV measurements taken per night [e.g., -nn 1].")
    parser.add_argument("-mhdn", "--minimum_hours_diff_night", type = float, default = None, help = "Minimum time in hours between consecutive observations happening in the same night [e.g., -mhdn 2.5].")

    # ONLY FOR Monotonic
    parser.add_argument("-cad", "--cadence", type = float, default = None, help = "Median cadence to observe the target in days [e.g., -cad 5].")

    # ONLY FOR Blocks
    parser.add_argument("-nb", "--number_blocks", type = int, default = None, help = "Maximum number of blocks [e.g., -nb 10]. Note that the simulated blocks can be below this number due to wheather loss and visibility conditions (if those are considered).")
    parser.add_argument("-db", "--duration_block", type = float, default = None, help = "Duration of a block in days [e.g., -db 10]. Number of consecutive nights that the target can be monitored during each block.")
    parser.add_argument("-gb", "--gap_blocks", type = float, default = None, help = "Gap between blocks in days [e.g., -gb 15]. This gap is from the last observation until the first one of a new block. If not given, it will be randomly chosen.")
    
    # ONLY FOR given JDs
    parser.add_argument("-jd", "--path_jd", type = str, default = None, help = "Path to the ASCII/CSV/TXT file with the BJDs to observe the target [e.g. -jd '/Desktop/jds_my_target.txt']. Format: column name must be 'jd' (the file can have more columns, but they will be ignored). Columns separated by commas (',').")


    #--------
    # Outputs
    #--------
    parser.add_argument("-dir", "--output_dir", default = None, help = "Output directory path [e.g., '/Desktop/moca_results/']")
    args = parser.parse_args()
    #--------------------------------------------


    # Initialise a dictionary to store all the inputs #
    dic_inputs = {}
    

    #-------------------#
    # Getting Star info #
    #-------------------#
    star = args.star_name
    ra = args.right_ascension
    dec = args.declination
    ms = args.stellar_mass
    st = args.spectral_type
    act_level = args.activity_level
    Prot = args.P_rotation
    Pmag = args.P_magnetic_cycle

    print('')
    print(' ')
    cprint("Star", 'green', attrs = ["underline"])
    if star == 'Unknown':
        text1 = colored('Name:', 'white')
        text2 = colored('Unknown', 'yellow')
        print(text1, text2)
        if ra and dec:
            print(f'RA = {ra} deg')
            print(f'DEC = {dec} deg')
            vis_bool_tar = True
        else:
            if not ra and dec:
                which_coord = 'RA'
            elif ra and not dec:
                which_coord = 'DEC'
            elif not ra and not dec:
                which_coord = 'RA, DEC'
            vis_bool_tar = False
            cprint(f'WARNING: {which_coord} and star name not provided', 'red')
    else:
        print(f'Name: {star}')
        vis_bool_tar = True
        if (not ra and dec) or (ra and not dec):
            if not ra:
                which_coord = 'RA'
            elif not dec:
                which_coord = 'DEC'
            cprint(f'WARNING: {which_coord} not provided. The coordinates are resolved by SIMBAD', 'red')
        if not (ra and dec):
            tar = Simbad.query_object(star)
            RA_tar, DEC_tar = tar['RA'][0], tar['DEC'][0]
            coord = SkyCoord(f'{RA_tar} {DEC_tar}', unit = (u.deg, u.deg))
            ra = coord.ra.value
            dec = coord.dec.value
        print(f'RA = {ra} deg')
        print(f'DEC = {dec} deg')


    if ms:
        print(f'Ms = {np.round(ms, 2)} Msun')
        if not st:
            st = st_to_mass.get_st(ms)
            text1 = colored(f'{st} V -', 'white')
            text2 = colored('inferred from Ms assuming main-sequence star', 'cyan')
            print(text1, text2)
        else:
            print(st, 'V')
        
        # check the given value for the stellar mass is reasonable
        if ms < 0.077:
            cprint('WARNING: The provided stellar mass is below 0.077 MSun, which is a substellar object!', 'red')
            cprint('Make sure -ms is in sollar mass units', 'red')
        elif ms > 2.5:
            cprint('WARNING: The provided Stellar Mass is above 2.5 MSun, which is hoter than a B9 ST!', 'red')
            cprint('Make sure -ms is in sollar mass units', 'red')

    elif st:
        ms = st_to_mass.get_ms(st) # ST is used to infer the stellar mass [Msun]
        print(st, 'V')
        text1 = colored(f'Ms = {np.round(ms, 2)} Msun -', 'white')
        text2 = colored('inferred from the ST assuming main-sequence star', 'cyan')
        print(text1, text2)
    else:
        text1 = colored('K7 V (0.646 Msun) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
        text1 = colored('NOTE:', 'blue')
        text2 = colored('you can provide either the stellar mass (-ms) or the spectral type (-st) of your target', 'white')
        print(text1, text2)
        ms = 0.646
        st = 'K7'

    if not act_level and act_level != 0:
        act_level = 0
        text1 = colored('Inactive star (0 m/s) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
    else:
        dic_act_level = {0: 'Inactive (0 m/s)', 1: 'Quiet (< 1 m/s)', 2: 'Moderate (1 - 5 m/s)', 3: 'Loud (5 - 20 m/s)'} # TD: modificate depending on the ST 
        if act_level not in dic_act_level.keys():
            cprint('ERROR: the activity level must be in the range of 0-3 from inactive to loud', 'red', attrs = ["bold"])
            sys.exit()
        print(f'Activity: {dic_act_level[act_level]}')
    
    if act_level != 0:
        if not Prot:
            Prot = 10 # TD: GENERAR RANDOM DEPENDIENDO DEL ST (SUAREZ-MASCAREÑO 2016)
            text1 = colored(f'Prot = {np.round(Prot, 2)} d -', 'white')
            text2 = colored('random based on the ST', 'yellow')
            print(text1, text2)
        else:
            print(f'Prot = {np.round(Prot, 2)} d')
        if not Pmag:
            Pmag = 150 # TD: GENERAR RANDOM DEPENDIENDO DEL ST (search)
            text1 = colored(f'Pmag = {np.round(Pmag, 2)} d -', 'white')
            text2 = colored('random based on the ST', 'yellow')
            print(text1, text2)
        else:
            print(f'Pmag = {np.round(Pmag, 2)} d')

    dic_inputs['star'] = star
    dic_inputs['ra'] = ra
    dic_inputs['dec'] = dec
    dic_inputs['ms'] = ms
    dic_inputs['st'] = st
    dic_inputs['act'] = act_level
    dic_inputs['Prot'] = Prot
    dic_inputs['Pmag'] = Pmag


    #------------------------#
    # Getting Planetary info #
    #------------------------#
    num_pl = args.number_planets
    per = args.orbital_period
    t0 = args.transiting_time
    ecc = args.eccentricity
    w = args.arg_periastron
    inc = args.inclination
    mp = args.planetary_mass
    rp = args.planetary_radius
    erp = args.unc_planetary_radius
    unit = args.units
    
    print(' ')
    print(' ')
    if not num_pl:
        cprint('ERROR: NUMBER OF PLANETS must be given (e.g., -np 1)', 'red', attrs = ["bold"])
        sys.exit()
    num_pl = num_pl
    if num_pl == 1:
        cprint("One planet", "green", attrs = ["underline"])
    else:
        cprint(f"{num_pl} planets", "green", attrs = ["underline"])


    if not unit:
        unit = 'e'
        text1 = colored('NOTE:', 'blue')
        text2 = colored('Earth units for the planetary mass/radius are used as default. You can change to Jupiter units with "-u j"', 'white')
        print(text1, text2)
    if unit == 'e':
        unit = 'Earth'
    elif unit == 'j':
        unit = 'Jupiter'
    else:
        cprint('ERROR: the unit for the planetary mass/radius must be "e" for Earth or "j" for Jupiter', 'red', attrs = ["bold"])
        sys.exit()


    if not per:
        per = np.random.uniform(0.3, 150., num_pl)
        persim = np.full(num_pl, True)
    elif 'nan' in per:
        persim = np.array([True if p == 'nan' else False for p in per])
        per = [np.round(np.random.uniform(0.3, 50.), 3) if p == 'nan' else p for p in per]
    else:
        persim = np.full(num_pl, False)
    per = list(map(float, per))
    if np.any(np.array(per) < 0):
        cprint('ERROR: The orbital period must be above 0 d', 'red', attrs = ["bold"])
        sys.exit()


    if not t0:
        t0 = np.random.uniform(2458849.5, 2460310.5, size = num_pl)
        t0sim = np.full(num_pl, True)
    elif 'nan' in t0:
        t0sim = np.array([True if t0i == 'nan' else False for t0i in t0])
        t0 = [np.round(np.random.uniform(2458849.5, 2460310.5), 8) if t0i == 'nan' else t0i for t0i in t0]
    else:
        t0sim = np.full(num_pl, False)
    t0 = list(map(float, t0))
    if np.any(np.array(t0) < 2378496):
        cprint('WARNING: Transiting time is before year 1800!', 'red')
    if np.any(np.array(t0) > Time(Time.now(), scale = 'utc').jd + 20*365): 
        cprint('WARNING: Transiting time is more than 20 years in the future!', 'red')
        

    if not ecc:
        ecc = truncnorm.rvs(loc = 0., scale = 0.2, a = 0, b = 1, size = num_pl)
        eccsim = np.full(num_pl, True)
    elif 'nan' in ecc:
        eccsim = np.array([True if e == 'nan' else False for e in ecc])
        ecc = [np.round(truncnorm.rvs(loc = 0., scale = 0.2, a = 0, b = 1), 3) if e == 'nan' else e for e in ecc]
    else:
        eccsim = np.full(num_pl, False)
    ecc = list(map(float, ecc))
    if np.any(np.array(ecc) < 0) or np.any(np.array(ecc) > 1):
        cprint('ERROR: Eccentricity must be between 0 and 1', 'red', attrs = ["bold"])
        sys.exit()


    if not w:
        w = np.random.uniform(0., 360., size = num_pl)
        wsim = np.full(num_pl, True)
    elif 'nan' in w:
        wsim = np.array([True if wi == 'nan' else False for wi in w])
        w = [np.round(np.random.uniform(0., 360.), 3) if wi == 'nan' else wi for wi in w]
    else:
        wsim = np.full(num_pl, False)
    w = list(map(float, w))
    if np.any(np.array(w) < 0) or np.any(np.array(w) > 360):
        cprint('ERROR: Argument of periastron must be between 0 and 360 deg', 'red', attrs = ["bold"])
        sys.exit()


    if not inc:
        inc = truncnorm.rvs(loc = 90., scale = 10., a = -90, b = 90, size = num_pl)
        incsim = np.full(num_pl, True)
    elif 'nan' in inc:
        incsim = np.array([True if i == 'nan' else False for i in inc])
        inc = [np.round(truncnorm.rvs(loc = 90., scale = 10., a = -90, b = 90), 3) if i == 'nan' else i for i in inc]
    else:
        incsim = np.full(num_pl, False)
    inc = list(map(float, inc))
    if np.any(np.array(inc) < 0) or np.any(np.array(inc) > 180):
        cprint('ERROR: Inclination must be between 0 and 180 deg', 'red', attrs = ["bold"])
        sys.exit()


    # Planetary mass
    mp_from_rp = np.array([])
    mpsim = np.array([])
    erpsim = np.array([])
    if not mp:
        mp = ['nan'] * num_pl
    if not rp:
        rp = ['nan'] * num_pl
    if not erp:
        erp = ['nan'] * num_pl

    for pl_i in range(num_pl):
        if mp[pl_i] == 'nan':
            if rp[pl_i] == 'nan': # No mass or radius is given for the planet
                mp_from_rp += False
                mpsim += True
                erpsim += False
                mp[pl_i] = np.random.uniform(0.8, 200.) 
            else: # estimate the planetary mass from the provided planetary radius using forecaster code (Chen & Kipping, 2017)
                mp_from_rp += True
                mpsim += False
                rp[pl_i] = float(rp[pl_i])
                if rp[pl_i] < 0:
                    cprint(f'ERROR: Planetary radius must be above 0 R{unit}', 'red', attrs = ["bold"])
                    sys.exit()
                if erp[pl_i] == 'nan':
                    erpsim += True
                    erp[pl_i] = 0.1 * rp[pl_i] # if no erp is given, it is used 10% Rp as default 
                else:
                    erp[pl_i] = float(erp[pl_i])
                    erpsim += False
                try:
                    mp[pl_i] = mr.Rstat2M(rp[pl_i], erp[pl_i], unit = unit)[0]
                except:
                    cprint(f'ERROR: Mass could not be estimated from the given Rp due to the eRp value associated (if not given, the default is 0.1 x Rp). Please, revise the value.', 'red', attrs = ["bold"])
                    sys.exit()
        else: # the planetary mass has been provided
            mp_from_rp = np.append(mp_from_rp, False)
            mpsim = np.append(mpsim, False)
            erpsim = np.append(erpsim, False)
        
        mp = list(map(float, mp))
        if np.any(np.array(mp) < 0):
            cprint(f'ERROR: Planetary mass must be above 0 M{unit}', 'red', attrs = ["bold"])
            sys.exit()


    # check that all parametes CALLED in the inputs have the same length 
    if not all([len(param) == num_pl for param in [per, t0, ecc, w, inc, mp]]):
        cprint('ERROR: not all the given planetary parameters have the same length', 'red', attrs = ["bold"])
        cprint('If you want to use a random value for a given parameter of a given planet use "nan" (e.g., -p 2.3 nan -e nan 0.3)', 'red')
        cprint('If a parameter is not called in the inputs, random values are generated for all the planets (e.g., "-np 2" will generate random values for all the orbital parameters for two planets)', 'red')
        sys.exit()

    for i_pl in range(num_pl):
        print_param_planet(i_pl, num_pl, per, persim, 'P', 'd', mp_from_rp)
        print_param_planet(i_pl, num_pl, t0, t0sim, 'T0', 'd', mp_from_rp)
        print_param_planet(i_pl, num_pl, ecc, eccsim, 'e', '', mp_from_rp)
        print_param_planet(i_pl, num_pl, inc, incsim, 'i', 'deg', mp_from_rp)
        print_param_planet(i_pl, num_pl, w, wsim, 'w', 'deg', mp_from_rp)
        #if rp[i_pl] != 'nan':
        #    print_param_planet(i_pl, num_pl, rp, mpsim, 'Rp', f'R{unit}', mp_from_rp, erp = erp, erpsim = erpsim)
        #print_param_planet(i_pl, num_pl, mp, mpsim, 'Mp', f'M{unit}', mp_from_rp)
        print(' ')


    dic_inputs['num_pl'] = num_pl
    dic_inputs['per'] = per
    dic_inputs['t0'] = t0
    dic_inputs['ecc'] = ecc
    dic_inputs['w'] = w
    dic_inputs['inc'] = inc
    dic_inputs['mp'] = mp
    dic_inputs['rp'] = rp
    dic_inputs['unit'] = unit

    #--------------------------------#
    # Getting the observational info #
    #--------------------------------#
    erv = args.rv_uncertainty
    obs = args.observatory
    wl = bool_str(args.weather_loss, '-wl')
    cwlp = args.custom_weather_loss_percentage
    id = args.init_date
    ed = args.end_date
    ma = args.minimum_altitude


    print('')
    cprint("Observational info", 'green', attrs = ["underline"])
    if not erv:
        erv = 3.0
        text1 = colored('Mean eRV: 3.0 m/s -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
    else:
        print(f'Mean eRV: {erv} m/s')

    
    if not obs:
        vis_bool_obs = False
        obs = 'other'
        obs_coords = None
        cprint('WARNING: No observatory provided', 'red')
    else: # TD: ask for weather loss statistics and add here the info
        vis_bool_obs = True
        dic_obs_coords_wl = {'CAHA': [[37.2236, -2.5461, 2168], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Large Binocular': [[32.7014, -109.8892, 3221], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'La Silla': [[-29.2609, -70.7315, 2400], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Lick': [[37.3411, -121.6431, 1283], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Mauna Kea': [[19.8230, -155.4694, 4205], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'McDonald': [[30.6717, -104.0217, 2075], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'OHP': [[43.9308, 5.7133, 650], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Paranal': [[-24.6276, -70.4051, 2635], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Roque Muchachos': [[28.7561, -17.8917, 2396], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]],
                            'Whipple': [[31.6811, -110.8781, 2340], [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]]}
        if obs not in dic_obs_coords_wl.keys():
            cprint(f'ERROR: {obs} is not considered as an observatory option. moca considers: {dic_obs_coords_wl.keys()}. If you detect an observatory hosting an echelle spectrograph missing in this list, please, contact us through GitHub or by email (obalsalobre@cab.inta-csic.es) to add it')
            sys.exit()
        else:
            print(f'Observatory: {obs}')
            obs_coords = dic_obs_coords_wl[obs][0]


    if wl != False and wl != True:
        text1 = colored('No weather loss considered -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
        wl = False
        wlp = None
    else:
        print(wl)
        if wl:
            if cwlp:
                print(f'An homogenous weather loss of the {cwlp}% is considered')
                wlp = [cwlp]*12
            elif obs == 'other':
                print(f'Weather loss is considered (homogenous loss of 30% over the year since no observatory is provided)')
                wlp = [30]*12
            else:
                print(f'Weather loss is considered for the {obs} observatory')
                wlp = dic_obs_coords_wl[obs][1]
        else:
            print('No weather loss considered')
            wlp = None
    

    if not id and not ed:
        id_jd = Time(Time.now(), scale = 'utc').jd + 1
        ed_jd = id_jd + 6 * 30
        id = Time(id_jd, format = 'jd').isot[:10]
        ed = Time(ed_jd, format = 'jd').isot[:10]
        text1 = colored(f'Init date: {id} (tomorrow) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
        text1 = colored(f'End date: {ed} (180 days apart from init date) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
    elif id and not ed:
        id_jd = Time(id, format = 'isot').jd
        ed_jd = id_jd + 6 * 30
        print(f'Init date: {id}')
        text1 = colored(f'End date: {ed} (180 days apart from init date) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
    elif not id and ed:
        ed_jd = Time(ed, format = 'isot').jd
        id_jd = ed_jd - 6 * 30
        print(f'End date: {ed}')
        text1 = colored(f'Init date: {id} (180 days before end date) -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
    else:
        id_jd = Time(id, format = 'isot').jd
        ed_jd = Time(ed, format = 'isot').jd
        print(f'Init date: {id}')
        print(f'End date: {ed}')


    if not vis_bool_tar and not vis_bool_obs:
        cprint('WARNING: Full visibility of the target is considered since the star (name or coordinates) and the observatory are not provided', 'red')
        vis = 'full'
    elif not vis_bool_tar and vis_bool_obs:
        cprint('WARNING: Full visibility of the target is considered since the star (name or coordinates) is not provided', 'red')
        vis = 'full'
    elif vis_bool_tar and not vis_bool_obs:
        cprint('WARNING: Full visibility of the target is considered since the observatory is not provided', 'red')
        vis = 'full'
    else:
        vis = 'check'
        if not ma:
            ma = 30
            text1 = colored('Minimum target altitude: 30 deg -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            print(f'Minimum target altitude: {ma} deg')


    dic_inputs['erv'] = erv
    dic_inputs['obs'] = obs
    dic_inputs['obs_coords'] = obs_coords
    dic_inputs['wl'] = wl
    dic_inputs['wlp'] = wlp
    dic_inputs['id'] = id
    dic_inputs['ed'] = ed
    dic_inputs['ma'] = ma
    dic_inputs['vis'] = vis


    #---------------------------#
    # Getting the strategy info #
    #---------------------------#
    sn = args.strategy_name
    nm = args.number_max_rvs
    nn = args.number_per_night_rvs
    mhdn = args.minimum_hours_diff_night
    cad = args.cadence
    nb = args.number_blocks
    db = args.duration_block
    gb = args.gap_blocks
    path_jd = args.path_jd
    
    
    print('')
    print(' ')
    cprint("Strategy", 'green', attrs = ["underline"])
    if not sn:
        text1 = colored('Name: monotonic -', 'white')
        text2 = colored('default', 'yellow')
        print(text1, text2)
        sn = 'mon'
        fsn = 'monotonic'
    else:
        dic_strategies = {'mon': 'monotonic', 'blo': 'blocks', 'jds': 'JDs Scheduled'}
        if sn not in dic_strategies.keys():
            cprint('ERROR: the given strategy name is not considered by moca. Options: ["mon", "blo", "jds"] for monotonic, blocks and JDs scheduled, respectively', 'red', attrs = ['bold'])
            sys.exit()
        fsn = dic_strategies[sn]
        print(f'Name: {fsn}')


    if sn == 'mon' or sn == 'blo':
        if not nm:
            nm = 30
            text1 = colored('Maximum Number of RVs = 30 -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if nm < 0:
                cprint('ERROR: Maximum Number of RVs must be above 0', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Maximum Number of RVs = {nm}')
        if not nn:
            nn = 1
            text1 = colored('Number of RVs per night = 1 -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if nn < 0:
                cprint('ERROR: Number of RVs per night must be above 0', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Number of RVs per night = {nn}')
            if not mhdn:
                mhdn = 30 / 60 # at least 30 minutes between observation in the same night
                text1 = colored(f'Minumum hours of difference in consecutive observations = {mhdn} -', 'white')
                text2 = colored('default', 'yellow')
                print(text1, text2)
            else:
                print(f'Minumum hours of difference in consecutive observations = {mhdn}')

    
    if sn == 'mon':
        if not cad:
            cad = 2
            text1 = colored('Cadence = 2 d -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if cad < 0:
                cprint('ERROR: Cadence must be above 0', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Cadence = {cad}')

    
    if sn == 'blo':
        if not nb:
            nb = 5
            text1 = colored('Number of blocks = 5 -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if nb < 0:
                cprint('ERROR: Number of blocks must be above 0', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Number of blocks = {nb}')
        if not db:
            db = 7
            text1 = colored('Duration per block = 7 d -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if db < 0:
                cprint('ERROR: Duration per block must be above 0', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Duration per block = {db} d')
        if not gb:
            gb = 10
            text1 = colored('Blocks gap 10 d -', 'white')
            text2 = colored('default', 'yellow')
            print(text1, text2)
        else:
            if gb < 1:
                cprint('ERROR: Blocks gap must be above 1', 'red', attrs = ['bold'])
                sys.exit()
            print(f'Blocks gap = {gb} d')

    
    if sn == 'jds':
        if not path_jd:
            cprint('ERROR: a path for the JDs scheduled must pe provided (e.g., -jd "Desktop/path_jd.txt"). It must be an ASCII/CSV/TXT file with a column of name "jd" collecting the observing times. Ifo more columns are provided, separation must be ","', 'red', attrs = ["bold"])
            sys.exit()
    else:
        jd_sch = None


    dic_inputs['sn'] = sn
    dic_inputs['fsn'] = fsn
    dic_inputs['nm'] = nm
    dic_inputs['nn'] = nn
    dic_inputs['mhdn'] = mhdn
    dic_inputs['cad'] = cad
    dic_inputs['nb'] = nb
    dic_inputs['db'] = db
    dic_inputs['gb'] = gb
    dic_inputs['path_jd'] = path_jd


    #--------------------------------#
    # Computing the RV semi-amplitude
    #--------------------------------#
    K = rv_functions.RV_semiamplitude(dic_inputs) # m/s
    dic_inputs['K'] = K


    #--------------------#
    # Getting output dir #
    #--------------------#
    directory = args.output_dir

    print(' ')
    print(' ')
    cprint("Outputs", "green", attrs = ["underline"])

    if not directory:
        directory = os.getcwd()
        cprint(f'Will be saved in Path: {directory} - default', 'yellow')
    else:
        print(f'Will be saved in Path: {directory}')
    dic_inputs['dir'] = directory
    print('')


    #--------------------------#
    # Save the simulation info #
    #--------------------------#
    dic_for_df = {k: [v] if isinstance(v, list) else v for k, v in dic_inputs.items()}
    df = pd.DataFrame(dic_for_df)
    tn = Time.now().isot[:16]

    full_dir = f'{directory}/moca_outputs/{star}'
    os.makedirs(full_dir, exist_ok = True)
    df.to_csv(full_dir+f'/sim_info_{tn}.csv', index = False, sep = ',')

    return dic_inputs