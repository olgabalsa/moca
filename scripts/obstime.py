import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_sun
from astropy.coordinates import SkyCoord
import astropy.units as u


def nights_lost_badweather(dic_inputs, jd):
    
    if not dic_inputs['wl']:
        return 0
    
    good_conditions = False
    n_lost_nights = 0
    while not good_conditions:
        likelihood_badw = np.random.uniform(0, 1)
        month = int(Time(jd, format = 'jd').isot[5:7])
        if likelihood_badw < dic_inputs['pwl'][month - 1]:
            good_conditions = False
            n_lost_nights += 1
        else:
            good_conditions = True
        
    return n_lost_nights


def lost_badweather(dic_inputs, jd):
    
    if not dic_inputs['wl']:
        return False
    
    likelihood_badw = np.random.uniform(0, 1)
    month = int(Time(jd, format = 'jd').isot[5:7])
    if likelihood_badw < dic_inputs['pwl'][month - 1]:
        return True
    return False


def range_tw(jd_arr, location, which_tw = 'both'):

    alt_astronomical_tw = -18 # deg

    range_tw_mo = []
    range_tw_ev = []
    for step, jdi in enumerate(jd_arr):
        time_format = Time(jdi, format = 'jd', scale = 'utc')
        sun_coord = get_sun(time_format).transform_to(AltAz(obstime = time_format, location = location)) # Alt, Az coordinates of the sun
        sun_alt = sun_coord.alt.value
        dif = sun_alt - alt_astronomical_tw
        if step == 0:
            prev_dif = dif
            continue # in the first step we define the variable 'prev_dif'
        if dif * prev_dif < 0: # check if there is a change in the sign of 'dif' in this step
            if dif > 0: # the altitude of the sun is increasing and therefore this is the morning twilight
                range_tw_mo = [jd_arr[step - 1], jdi]
            else: # the altitude of the sun is decreasing and therefore this is the evening twilight
                range_tw_ev = [jd_arr[step - 1], jdi]

        # check if the target twilight has been found to stop the loop
        if (len(range_tw_mo) != 0) and (len(range_tw_ev) != 0) and (which_tw == 'both'):
            break
        elif (len(range_tw_mo) != 0) and (which_tw == 'mo'):
            break
        elif (len(range_tw_ev) != 0) and (which_tw == 'ev'):
            break
        
        prev_dif = dif # redifine for the next step

    return range_tw_mo, range_tw_ev


def compute_tw(jd_date, obs_loc):
    
    # first approximation of the twilightls by using big steps of 1 h
    jd_init = jd_date
    jd_end = jd_init + 1.0
    jd_h_step = np.linspace(jd_init, jd_end, 25)
    range_tw_mo, range_tw_ev = range_tw(jd_h_step, obs_loc)

    # we fine tune the twilightls by using smaller steps of 1 min
    #print('MO',range_tw_mo, range_tw_mo)
    jd_tw_mo_array = np.linspace(range_tw_mo[0], range_tw_mo[1], 61) 
    small_range_tw_mo, _ = range_tw(jd_tw_mo_array, obs_loc, which_tw = 'mo')
    tw_mo = np.median(small_range_tw_mo)
    
    #print('EV',range_tw_ev, range_tw_ev)
    jd_tw_ev_array = np.linspace(range_tw_ev[0], range_tw_ev[1], 61) 
    _, small_range_tw_ev = range_tw(jd_tw_ev_array, obs_loc, which_tw = 'ev')
    tw_ev = np.median(small_range_tw_ev)

    return tw_mo, tw_ev


def find_jd_visible(jd, dic_inputs, repeat_night):

    if dic_inputs['vis'] == 'full':
        return jd + np.random.normal(0, 0.4)
    
    lat, lon, alt = dic_inputs['obs_coords']
    obs_loc = EarthLocation(lat = lat * u.deg, lon = lon * u.deg, height = alt * u.m)

    try:
        tw_ev, tw_mo = compute_tw(jd, obs_loc)
    except: # tw matching the jd so it explodes
        jd += 0.01
        tw_ev, tw_mo = compute_tw(jd, obs_loc)
    if tw_ev > tw_mo: # the jd corresponds with an already started night, so we compute the start of it 
        try:
            tw_ev, _ = compute_tw(jd - 1, obs_loc)
        except:
            jd += 0.01
            tw_ev, _ = compute_tw(jd - 1, obs_loc)
    if repeat_night:
        if not (jd > tw_ev and jd < tw_mo): # the night is finished
            return np.nan
        else:
            jd_init = jd
    else:
        jd_init = tw_ev
    
    coord = SkyCoord(f'{dic_inputs["ra"]} {dic_inputs["dec"]}', unit = (u.deg, u.deg))
    jd_night = np.linspace(jd_init, tw_mo, 30)
    star_alt = coord.transform_to(AltAz(location = obs_loc, obstime = Time(jd_night, format = 'jd'))).alt.value
    index_mask_visible = np.where(star_alt > dic_inputs['ma'])[0]
    if len(jd_night[index_mask_visible]) == 0:
        return np.nan
    return jd_night[index_mask_visible][0] # select the first chance to observe the target


def find_new_jd(dic_inputs, jd_i, measure_of_the_night, nights_current_block):

    jd_next_visible = np.nan # initialize variable
    sn = dic_inputs['sn']
    
    while np.isnan(jd_next_visible):
        #---- Next date based on the strategy ----#
        if sn == 'jds':
            arg_last_jd = np.where(dic_inputs['jd_sch'] == jd_i)[0][0]
            try:
                expected_jd = dic_inputs['jd_sch'][arg_last_jd + 1]
                n_lost_nights = nights_lost_badweather(dic_inputs, expected_jd)
                jd_next_visible = dic_inputs['jd_sch'][arg_last_jd + 1 + n_lost_nights]
            except:
                return np.nan, measure_of_the_night, nights_current_block
        
        else:

            # Observe again during the same night?
            lost_night = lost_badweather(dic_inputs, jd_i)
            if not lost_night and dic_inputs['nn'] > measure_of_the_night:
                repeat_night = True
                measure_of_the_night += 1
                jd_next = jd_i + np.random.uniform(dic_inputs['mhdn']/24 + (1.5+dic_inputs['mhdn'])/24)
            else:
                measure_of_the_night = 1
                repeat_night = False

            if sn == 'mon' and not repeat_night:
                lost_night = True # init
                n_lost_nights = 0
                while lost_night:
                    jd_next = jd_i + (n_lost_nights + 1) * dic_inputs['cad']
                    lost_night = lost_badweather(dic_inputs, jd_next)
                    n_lost_nights += 1

            if sn == 'blo' and not repeat_night:
                lost_night = True # init
                n_lost_nights = 0
                while lost_night:
                    jd_next = jd_i + 1
                    nights_current_block += 1
                    if nights_current_block > dic_inputs['db']:
                        jd_next = jd_i + dic_inputs['gb']
                        nights_current_block = 1
                    lost_night = lost_badweather(dic_inputs, jd_i)
                    n_lost_nights += 1

            jd_next_visible = find_jd_visible(jd_next, dic_inputs, repeat_night)
            if jd_next_visible > Time(dic_inputs['ed'], format = 'isot', scale = 'utc').jd:
                return np.nan, measure_of_the_night, nights_current_block
            jd_i = jd_next # to re-do the search in case it is not visible

    return jd_next_visible, measure_of_the_night, nights_current_block