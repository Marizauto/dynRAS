import numpy as np


def Weight(Tank, params, t):
    dweight =  3 * params.TGC* params.T /((24*60*60)*1000) * (
            params.w0**(1 / 3) + params.TGC * params.T * t / ((24*60*60) * 1000)) ** 2
    dweight = dweight

    dBiomass = (3 * params.TGC* params.T /((24*60*60)*1000) * (
            params.w0**(1 / 3) + params.TGC * params.T * t / ((24*60*60) * 1000)) ** 2) * Tank.Fish_number
    dY = np.array([dweight, dBiomass])

    return dY