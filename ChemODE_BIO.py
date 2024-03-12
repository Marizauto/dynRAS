import numpy as np
from Growth_Bacteria import growth_bact
def chemODE_BIO(Tank,params,t):
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    day = t / 60 / 60 / 24
#simulation case 1
    if day <= 14:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day > 14 and day <= 28:
        alkalinity = 100/50.04
        pH=-np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 28 and day <= 42:
        alkalinity = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 42 and day <= 56:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 56 and day <= 70:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day > 70 and day <= 84:
        alkalinity = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 84 and day <= 98:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 98 and day <= 112:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0
    elif day > 112 and day <= 126:

        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 126 and day <= 140:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 140 and day <= 154:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0

    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = params.k1*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB*1/(params.YAOB/((1/14)/(7/61.1)))/Tank.V+add_HCO3+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = params.k3-params.k_3*Tank.H*Tank.OH+add_OH+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =(-params.k4*Tank.NH4)+params.k_4*Tank.NH3*Tank.H+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*(1/params.YAOB)+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=+params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*1/params.YAOB-params.muNOB*(Tank.NO2/(Tank.NO2+params.KNO2))*Tank.NOB/Tank.V*1/params.YNOB+(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dAOB=g_AOB-d_AOB
    dNOB=g_NOB-d_NOB
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2,dAOB,dNOB])
    return dY


def chemODE_BIO_HCO3(Tank,params,t):
    #simulation case 2
    # HCO3 only
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    day = t / 60 / 60 / 24
    if day < 14:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 14 and day < 28:
        desired_HCO3_concentration = 100/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 28 and day < 42:
        desired_HCO3_concentration = 100/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 42 and day < 56:
        desired_HCO3_concentration = 70/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 56 and day < 70:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 70 and day < 84:
        desired_HCO3_concentration =100/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 84 and day < 98:
        desired_HCO3_concentration = 70/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day >= 98 and day < 112:
        desired_HCO3_concentration = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0
    elif day >= 112 and day < 126:
        desired_HCO3_concentration = 70/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0
    elif day >= 126 and day <= 140:
        desired_HCO3_concentration = 70/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_HCO3_concentration) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0
    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = params.k1*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB*1/(params.YAOB/((1/14)/(7/61.1)))/Tank.V+add_HCO3+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = params.k3-params.k_3*Tank.H*Tank.OH+add_OH+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =(-params.k4*Tank.NH4)+params.k_4*Tank.NH3*Tank.H+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*(1/params.YAOB)+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=+params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*1/params.YAOB-params.muNOB*(Tank.NO2/(Tank.NO2+params.KNO2))*Tank.NOB/Tank.V*1/params.YNOB+(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dAOB=g_AOB-d_AOB
    dNOB=g_NOB-d_NOB
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2,dAOB,dNOB])
    return dY


def chemODE_BIO_NaOH(Tank, params, t):
    #simulation case 2
    # NaOH only
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    day = t / 60 / 60 / 24
    if day < 14:
        alkalinity = 200/50.04
        pH=-np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 14 and day < 28:
        alkalinity = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0

    elif day >= 28 and day < 42:
        alkalinity = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 42 and day < 56:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0


    elif day >= 56 and day < 70:
        alkalinity = 200/50.04
        pH=-np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        if desired_concentration_min < 10**-8*10**3:
            desired_concentration_min = 10**-8*10**3
            feedback_term_OH = (1 + (desired_concentration_min / Tank.H)**10)
            add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        else:
            desired_concentration_min = desired_concentration_min
            feedback_term_OH = (1 + (desired_concentration_min / Tank.H)**10)
            add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0

    elif day >= 70 and day < 84:
        alkalinity = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H)**10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 84 and day < 98:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 98 and day < 112:
        alkalinity =200/50.04
        pH=-np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        if desired_concentration_min < 10**-8*10**3:
            desired_concentration_min = 10**-8*10**3
            feedback_term_OH = (1 + (desired_concentration_min / Tank.H)**10)
            add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        else:
            desired_concentration_min = desired_concentration_min
            feedback_term_OH = (1 + (desired_concentration_min / Tank.H)**10)
            add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 112 and day < 126:

        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day >= 126 and day <= 140:
        alkalinity = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((alkalinity/(Tank.CO2aq)))
        desired_concentration_min=(10**(-pH))*10**3
        feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = params.k1*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB*1/(params.YAOB/((1/14)/(7/61.1)))/Tank.V+add_HCO3+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = params.k3-params.k_3*Tank.H*Tank.OH+add_OH+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =(-params.k4*Tank.NH4)+params.k_4*Tank.NH3*Tank.H+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*(1/params.YAOB)+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=+params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*1/params.YAOB-params.muNOB*(Tank.NO2/(Tank.NO2+params.KNO2))*Tank.NOB/Tank.V*1/params.YNOB+(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dAOB=g_AOB-d_AOB
    dNOB=g_NOB-d_NOB
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2,dAOB,dNOB])
    return dY
def chemODE_BIO_pH_control(Tank, params, t):
    #set the pH with either NaOH or HCO3
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    add_OH = 0 #comment or uncomment this line to turn on or off pH control using NaOH
    desired_pH = 7
    #desired_concentration_min = (10 ** -desired_pH) * 10 ** 3 #dosing with OH
    #feedback_term_OH = (1 + (desired_concentration_min / Tank.H) ** 10) #dosing with OH
    #add_OH = (0.0001* Tank.dosing_OH.OH / Tank.V) / feedback_term_OH #dosing with OH
    desired_concentration_min_HCO3= (10**(desired_pH-6.1))*Tank.CO2aq #dosing with HCO3
    feedback_term_HCO3 = (1 + (Tank.HCO3 / desired_concentration_min_HCO3)**4) #dosing with HCO3
    add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V)/feedback_term_HCO3 #dosing with HCO3
    # add_HCO3 = 0 #comment or uncomment this line to turn on or off pH control using HCO3
    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    dHCO3 = params.k1*Tank.CO2aq-(params.k_1*Tank.H)*Tank.HCO3+params.k_H*Tank.H*Tank.CO32-params.kH*Tank.HCO3-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB*1/(params.YAOB/((1/14)/(7/61.1)))/Tank.V+add_HCO3+(Tank.Tminus.HCO3-Tank.HCO3)*params.F/Tank.V
    dCO32 =params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+(Tank.Tminus.CO32-Tank.CO32)*params.F/Tank.V
    dH = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H+params.k1*Tank.CO2aq-params.k_1*Tank.H*Tank.HCO3+params.kH*Tank.HCO3-params.k_H*Tank.H*Tank.CO32+params.k3-params.k_3*Tank.H*Tank.OH+(Tank.Tminus.H-Tank.H)*params.F/Tank.V
    dOH = params.k3-params.k_3*Tank.H*Tank.OH+add_OH+(Tank.Tminus.OH-Tank.OH)*params.F/Tank.V
    dNH4 =(-params.k4*Tank.NH4)+params.k_4*Tank.NH3*Tank.H+(Tank.Tminus.NH4-Tank.NH4)*params.F/Tank.V
    dNH3 = params.k4*Tank.NH4-params.k_4*Tank.NH3*Tank.H-params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*(1/params.YAOB)+(Tank.Tminus.NH3-Tank.NH3)*params.F/Tank.V
    dNO2=+params.muAOB*(Tank.NH3/(params.KNH3+Tank.NH3))*(Tank.HCO3/(params.KAlk+Tank.HCO3))*Tank.AOB/Tank.V*1/params.YAOB-params.muNOB*(Tank.NO2/(Tank.NO2+params.KNO2))*Tank.NOB/Tank.V*1/params.YNOB+(Tank.Tminus.NO2-Tank.NO2)*params.F/Tank.V
    dAOB=g_AOB-d_AOB
    dNOB=g_NOB-d_NOB
    dY = np.array([dCO2aq, dHCO3, dCO32, dH, dOH, dNH4, dNH3,dNO2,dAOB,dNOB])
    return dY