import numpy as np
import math
from TAN_prod import NH4_rate
def chemODE_FISH(Tank, params,t):
    day=t/60/60/24
    print(day)
    interval=24*(day-math.floor(day))
    if interval<12:
        CO2P =((((62.5*(Tank.Fish_weight*10**-3)**-0.3)*1.06**14*(44/32)/(60*60))*Tank.Fish_Biomass/1000)/44/Tank.V)*1.037**interval
    else:
        CO2P = ((((62.5*(Tank.Fish_weight*10**-3)**-0.3)*1.06**14*(44/32)/(60*60))*Tank.Fish_Biomass/1000)/44/Tank.V)* 1.037** (24-interval)
    TAN_total=NH4_rate(Tank,t)
    K = 2.73*10**-10

    # Calculate [H+] from pH
    H_concentration = Tank.H*10**-3

    # Calculate the fraction of NH3 in the total TAN
    NH3_fraction = 1 / (1 + (H_concentration/ K))

    # The fraction of NH4+ will be the remainder
    NH4_fraction = 1 - NH3_fraction

    # Calculate the amounts of NH3 and NH4+
    NH3_total = NH3_fraction * TAN_total
    NH4_total = NH4_fraction * TAN_total

    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+CO2P+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = (params.k1)*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = +params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = +params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =-params.k4*Tank.NH4+params.k_4*Tank.NH3*Tank.H+NH4_total/(24*60*60)/Tank.V+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = +params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+NH3_total/(24*60*60)/Tank.V+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2])

    return dY