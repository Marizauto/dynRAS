import numpy as np
from Growth_Bacteria import growth_bact
def chemODE_BIO(Tank,params,t):
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    day = t / 60 / 60 / 24
    #K1=(johnson 1982)
# simulation Jafari et al. 2024
    if day <= 14:
        T_HCO3 = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day > 14 and day <= 28:
        T_HCO3 = 100/50.04
        pH=-np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 28 and day <= 42:
        T_HCO3 = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 42 and day <= 56:
        T_HCO3 = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 56 and day <= 70:
        T_HCO3 = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH= 0
    elif day > 70 and day <= 84:
        T_HCO3 = 100/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 84 and day <= 98:
        T_HCO3 = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 98 and day <= 112:
        T_HCO3 = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
        add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
        add_OH = 0
    elif day > 112 and day <= 126:

        T_HCO3 = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 126 and day <= 140:
        T_HCO3 = 70/50.04
        pH= -np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
        T_H=(10**(-pH))*10**3
        feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
        add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
        add_HCO3 = 0
    elif day > 140 and day <= 154:
        T_HCO3 = 200/50.04
        feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
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
    #simulation case 1
    # HCO3 only
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    T_HCO3= 200/50.04
    #desired_HCO3_concentration = 70/ 50.04 uncomment to change the desired concentration
    feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
    add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3
    add_OH= 0
    Tank.dosing_amount_HCO3.append(add_HCO3)
    Tank.dosing_time.append(t)
    Tank.dosing_amount_OH.append(add_OH)

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
    T_HCO3 = 200/50.04
    #T_HCO3 = 70/ 50.04 uncomment to change the desired concentration
    pH=-np.log10(0.79*10**-6)+np.log10((T_HCO3/(Tank.CO2aq)))
    T_H=(10**(-pH))*10**3
    if T_H<10**-8*10**3:
        T_H=10**-8*10**3
    feedback_term_OH = (1 + (T_H / Tank.H) ** 2)
    add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / feedback_term_OH
    add_HCO3 = 0
    Tank.dosing_time.append(t)
    Tank.dosing_amount_OH.append(add_OH)
    Tank.dosing_amount_HCO3.append(add_HCO3)
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
   #add_OH = 0 #comment or uncomment this line to turn on or off pH control using NaOH
    desired_pH = 7
    T_H = (10 ** -desired_pH) * 10 ** 3 #dosing with OH
    add_OH = 0  # comment or uncomment this line to turn on or off pH control using NaOH
    feedback_term_OH = (1 + (T_H / Tank.H) ** 10) #dosing with OH
    #add_OH = (0.0001* Tank.dosing_OH.OH / Tank.V) / feedback_term_OH #dosing with OH
    Tank.dosing_time.append(t)
    Tank.dosing_amount_OH.append(add_OH)
    T_HCO3 = (10 ** (desired_pH - (-np.log10(0.79*10**-6)))) * Tank.CO2aq  # dosing with HCO3
    feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)  # dosing with HCO3
    add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / feedback_term_HCO3  # dosing with HCO3
    #add_HCO3 = 0 #comment or uncomment this line to turn on or off pH control using HCO3
    Tank.dosing_amount_HCO3.append(add_HCO3)
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
def chemODE_BIO_alk_control(Tank, params, t):
    #set alkalnity with either NaOH or HCO3 based on the CO2 in the system
    Growth=growth_bact(Tank,params)
    g_AOB=Growth[0]
    g_NOB=Growth[1]
    d_AOB=Growth[2]
    d_NOB=Growth[3]
    T_HCO3 = 200/ 50.04
    CO2_Threshold=10/44.1
    pH = -np.log10(0.79 * 10 ** -6) + np.log10((T_HCO3 / (Tank.CO2aq)))
    T_H = (10 ** (-pH)) * 10 ** 3
    feedback_term_OH = (1 + (T_H / Tank.H) ** 10)
    feedback_CO2_OH=(1+(CO2_Threshold/Tank.CO2aq)**10)
    feedback_CO2_HCO3=(1+(Tank.CO2aq/CO2_Threshold)**10)
    add_OH = (0.0001 * Tank.dosing_OH.OH / Tank.V) / (feedback_term_OH*feedback_CO2_OH)
    feedback_term_HCO3 = (1 + (Tank.HCO3 / T_HCO3) ** 10)
    add_HCO3 = (0.0001 * Tank.dosing_HCO3.HCO3 / Tank.V) / (feedback_term_HCO3*feedback_CO2_HCO3) # dosing with HCO3
    dCO2aq = (params.k_1*Tank.H)*Tank.HCO3-(params.k1)*Tank.CO2aq+(Tank.Tminus.CO2aq-Tank.CO2aq)*params.F/Tank.V
    Tank.dosing_time.append(t)
    Tank.dosing_amount_OH.append(add_OH)
    Tank.dosing_amount_HCO3.append(add_HCO3)
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