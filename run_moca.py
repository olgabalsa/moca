import scripts.inputs as inputs # type: ignore
import scripts.rv_functions as rv_functions # type: ignore
import scripts.outputs as outputs # type: ignore
import scripts.obstime as obstime # type: ignore
from termcolor import cprint
import numpy as np
import pandas as pd
from astropy.time import Time
import sys


dic_inputs = inputs.get()


# Initilizing the observing dates
if dic_inputs['sn'] == 'jds':
    jd_sch = pd.read_csv(dic_inputs['path_jd'], sep = ',').jd.values
    dic_inputs['jd_sch'] = jd_sch
    jd_i = jd_sch[0]
    if jd_i < Time(dic_inputs['id'], format = 'isot', scale = 'utc').jd:
        cprint('ERROR: the first day of your schedule is BEFORE the beginning of the observations. You should increase the -id or change your schedule).', 'red')
        sys.exit()
else:
    jd_i = Time(dic_inputs['id'], format = 'isot', scale = 'utc').jd + 0.5


# Computing the simulated dataset
mock_timeseries = {'jd': np.array([]), 'rv': np.array([]), 'erv': np.array([])}


# inits of some variables
measure_of_the_night = 1
nights_current_block = 1
number_lost_nights = 0


print('Mocking the RVs in process...')
for ndp in range(dic_inputs['nm']):

    rv_i = rv_functions.rv_model(np.array([jd_i]), dic_inputs)
    noise_i = np.random.normal(0., dic_inputs['erv']) # jitter cosidered to be of the same magnitude as noise
    rv_i += noise_i
    erv_i = np.random.normal(dic_inputs['erv'], 0.2 * dic_inputs['erv'])

    mock_timeseries['jd'] = np.append(mock_timeseries['jd'], jd_i)
    mock_timeseries['rv'] = np.append(mock_timeseries['rv'], rv_i)
    mock_timeseries['erv'] = np.append(mock_timeseries['erv'], erv_i)

    # Find the next JD to observe
    jd_i, measure_of_the_night, nights_current_block = obstime.find_new_jd(dic_inputs, jd_i, measure_of_the_night, nights_current_block)
    print(f'Mock RV #{ndp+1} \u2713')

    if np.isnan(jd_i):
        break


# Outputs
print('Saving the dataset')
outputs.save_ascii(mock_timeseries, dic_inputs)
print('Saving figures')
outputs.do_plot(mock_timeseries, dic_inputs)
print('Done!')
print('')