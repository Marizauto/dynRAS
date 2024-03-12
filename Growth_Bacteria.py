
def growth_bact(Tank,params):

#4.	Aerobic growth of AOB
    g_AOB=params.muAOB*((Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3)))*Tank.AOB*(Tank.B-Tank.AOB-Tank.NOB)
#5.	Aerobic growth of NOB
    g_NOB=params.muNOB*((Tank.NO2/(params.KNO2+Tank.NO2)))*Tank.NOB*(Tank.B-Tank.AOB-Tank.NOB)

#7. Decay of AOB
    d_AOB=params.rhoAOB*Tank.AOB
#	8. Decay of NOB
    d_NOB=params.rhoNOB*Tank.NOB

    growth= [g_AOB, g_NOB,d_AOB,d_NOB]

    return growth