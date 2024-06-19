#chemical ODEs for Degasser
#k=forward rate
#k_=reverse rate
import numpy as np
def chemODE_DGS(Tank, params,t):
    day = t / 60 / 60 / 24
    r=0.0015
    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq-Tank.CO2aq*r+(Tank.Tplus.CO2aq-Tank.CO2aq)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = (params.k1)*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3+(Tank.Tplus.HCO3-Tank.HCO3)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tplus.CO32-Tank.CO32)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tplus.H-Tank.H)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tplus.OH-Tank.OH)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =-params.k4*Tank.NH4+params.k_4*Tank.NH3*Tank.H+(Tank.Tplus.NH4-Tank.NH4)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+(Tank.Tplus.NH3-Tank.NH3)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=+(Tank.Tplus.NO2-Tank.NO2)*Tank.Tplus.exchange_rate/Tank.V+(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2])
    return dY