import scripts.inputs as inputs # type: ignore
import scripts.rv_functions as rv_functions # type: ignore
import scripts.outputs as outputs # type: ignore
import numpy as np
from astropy.time import Time
import sys


dic_inputs = inputs.get()


# Initilizing the observing dates
if dic_inputs['sn'] == 'jds':
    jd = dic_inputs['jd_sch']
else:
    jd_i = Time(Time.now(), scale = 'utc').jd


# Architecture of the planetary system
ecc = dic_inputs['ecc']
w = dic_inputs['w']
c = ecc * np.cos(np.radians(w))
d = ecc * np.sin(np.radians(w))










# Observational info


# Strategy info
nrv = dic_inputs['nm']
cad = dic_inputs['cad']

#########
# Star info
Prot = dic_inputs['Prot']
Pmag = dic_inputs['Pmag']


# Let's create those Mock RVs!!!
mocha_rvs = np.array([])
mocha_ervs = np.array([])
mocha_jds = np.array([])
for ndp in range(nrv):

    rv_i = rv_functions.rv_model(np.array([jd_i]), dic_inputs['num_pl'], dic_inputs['K'], dic_inputs['per'], dic_inputs['t0'], c, d)
    noise_i = np.random.normal(0., dic_inputs['erv']) # jitter cosidered to be of the same magnitude as noise
    rv_i += noise_i
    erv_i = np.random.normal(erv, 0.2 * dic_inputs['erv'])

    # store the mock values in arrays
    mocha_rvs = np.append(mocha_rvs, rv_i)
    mocha_ervs = np.append(mocha_ervs, erv_i)
    mocha_jds = np.append(mocha_jds, jd_i)

    jd_i += np.random.normal(cad, 0.15) # new date for the next iteration


# Outputs
dir = dic_inputs['dir']
star = dic_inputs['star']
outputs.save_ascii(mocha_rvs, mocha_ervs, mocha_jds, dir, star)
outputs.plot_RV_jd(mocha_jds, num_pl, mocha_rvs, mocha_ervs, K, per, t0, c, d, Prot, Pmag, dir, star)

print('Simulation finished!')
print('')