import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ChemODE_BIO import chemODE_BIO #simulation case 1 modify uncomment function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_HCO3 #simulation case 2 HCO3 dosing only uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_NaOH #simulation case 2 NaOH dosing only uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_BIO import chemODE_BIO_pH_control #set a pH threshold and dose NaOH or HCO3 to maintain the pH at the setpoint uncomment the function for the right hand side of for the ODE of the biofilter class to run it
from ChemODE_Fish import chemODE_FISH
from scipy.integrate import solve_ivp
from params import params
from Fish_growth import Weight
from Biomass_function import Biomass
from ChemODE_DGS import chemODE_DGS
class Sump:
    def __init__(self, CO2aq, HCO3, CO32, H, OH,NH4, NH3, NO2,V,exchange_rate):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.V = V
        self.NO2 = NO2
        self.exchange_rate = exchange_rate

class Dosing_pump:
    def __init__(self, OH, dosing_compartment,V):
        self.OH = OH
        self.dosing_compartment = dosing_compartment
        self.V=V

class Dosing_pump_HCO3:
    def __init__(self, HCO3, dosing_compartment,V):
        self.HCO3 = HCO3
        self.dosing_compartment = dosing_compartment
        self.V=V

class fish_tank:
    Fish_Biomass_Max = 50000

    def __init__(self, CO2aq, HCO3, CO32, H, OH, NH4, NH3,NO2, Tminus, V,
                 Fish_weight, Fish_number, Fish_Biomass):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.NO2 = NO2
        self.Tminus = Tminus
        self.V = V
        self.Fish_weight = Fish_weight
        self.Fish_number = Fish_number
        self.Fish_Biomass = Fish_Biomass

    def RHS(self, params, t):
        dY1 = chemODE_FISH(self, params,t)
        dY2 = Weight(self, params, t)
        dY = np.append(dY1, dY2)
        return dY




class biofilter:
    def __init__(self,CO2aq, HCO3, CO32, H, OH, NH4, NH3, NO2, Tminus,dosing_OH,dosing_HCO3, V, AOB,
                 NOB,B):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH =OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.NO2 = NO2
        self.Tminus = Tminus
        self.dosing_OH = dosing_OH
        self.dosing_HCO3 = dosing_HCO3
        self.V = V
        self.AOB = AOB
        self.NOB = NOB
        self.B=B

    def RHS(self, params, t):
        #dY = chemODE_BIO(self, params,t) #uncomment to run and comment the three other lines (simulation case 1)
        dY = chemODE_BIO_HCO3(self, params, t) #uncomment to run and comment the three other lines (simulation case 2_HCO3 dosing only)
        #dY = chemODE_BIO_NaOH(self, params, t) #uncomment to run and comment the three other lines (simulation case 2_NaOH dosing only)
        #dY=chemODE_BIO_pH_control(self,params, t) #uncomment to run and comment the three other lines set a threshold for the pH and control it using eiter OH or HCO3
        return dY


class degasser:
    def __init__(self, CO2aq, HCO3, CO32, H, OH, NH4, NH3, NO2, Tminus, V,Tplus):
        self.CO2aq = CO2aq
        self.HCO3 = HCO3
        self.CO32 = CO32
        self.H = H
        self.OH = OH
        self.NH4 = NH4
        self.NH3 = NH3
        self.Tminus = Tminus
        self.V = V
        self.NO2 = NO2
        self.Tplus=Tplus

    def RHS(self, params, t):
        dY1 = chemODE_DGS(self, params,t)
        return dY1


def chem(t, S, params, FishTank, B1, DGS):
    # Assigning values from the array S to the  FishTank
    FishTank.CO2aq = S[0]
    FishTank.HCO3 = S[1]
    FishTank.CO32 = S[2]
    FishTank.H = S[3]
    FishTank.OH = S[4]
    FishTank.NH4 = S[5]
    FishTank.NH3 = S[6]
    FishTank.NO2 = S[7]
    FishTank.Fish_weight = S[8]
    FishTank.Fish_Biomass = S[9]

    # Assigning values to B1
    B1.CO2aq = S[10]
    B1.HCO3 = S[11]
    B1.CO32 = S[12]
    B1.H = S[13]
    B1.OH = S[14]
    B1.NH4 = S[15]
    B1.NH3 = S[16]
    B1.NO2 = S[17]
    B1.AOB = S[18]
    B1.NOB = S[19]

    # Assigning values to DGS
    DGS.CO2aq = S[20]
    DGS.HCO3 = S[21]
    DGS.CO32 = S[22]
    DGS.H = S[23]
    DGS.OH = S[24]
    DGS.NH4 = S[25]
    DGS.NH3 = S[26]
    DGS.NO2 = S[27]

    # Compute RHS for each component
    dY1 = FishTank.RHS(params, t)
    dY2 = B1.RHS(params, t)
    dY3 = DGS.RHS(params, t)
    dY4 = np.append(dY1, dY2)
    dY = np.append(dY4, dY3)

    return dY


# initial condition Fish tank


CO2aq_FT = 10/44.1 #10 mg/L / 44.1 g/mol 10/44.1 mmol/L
HCO3_FT =200/61.01 #(mg/l)==> mmol.l-1
CO32_FT =10**-6*10**3
H_FT = 10**-7.6*10**3
OH_FT = 10**-6.4*10**3
NH4_FT =0
NH3_FT =0
NO2_FT =0

Fishweight_FT = params.w0

Fish_number_FT = 525

Fish_Biomass_FT = Fish_number_FT * Fishweight_FT

# initial condition B1


CO2aq_B1 =10/44.1
HCO3_B1 =200/61.01
CO32_B1 =10**-6*10**3
H_B1 =10**-7.6*10**3
OH_B1 = 10**-6.5*10**3
NH4_B1 =0
NH3_B1 =0
NO2_B1 =0
AOB_B1 = 960000
NOB_B1 = 480000
B=0.3*80000*10**3*0.05
# initial condition DGS
CO2aq_DGS =10/44.1
HCO3_DGS =200/61.01
CO32_DGS=10**-6*10**3
H_DGS =10**-7.6*10**3
OH_DGS =10**-6.4*10**3
NH4_DGS =0
NH3_DGS =0
NO2_DGS =0
NO3_DGS =40/62


# initial values for the dosing pump
#pump_OH
OH_dosing=2000 #mmol/l
#pump_HCO3
HCO3_dosing=2000
#value for the sump
CO2aq_Sump =0
HCO3_Sump =70/61.01
CO32_Sump=10**-6*10**3
H_Sump =10**-7.5*10**3
OH_Sump =10**-6.5*10**3
NH4_Sump =0
NH3_Sump =0
NO2_Sump =0
V_Sump = 50



FishTank = fish_tank( CO2aq_FT, HCO3_FT, CO32_FT, H_FT, OH_FT,NH4_FT, NH3_FT,
                     NO2_FT, [], 1000, Fishweight_FT, Fish_number_FT, Fish_Biomass_FT)
B1 = biofilter( CO2aq_B1, HCO3_B1, CO32_B1, H_B1, OH_B1, NH4_B1, NH3_B1,
             NO2_B1, FishTank,[],[], 800, AOB_B1, NOB_B1, B)
DGS = degasser( CO2aq_DGS, HCO3_DGS, CO32_DGS, H_DGS, OH_DGS, NH4_DGS, NH3_DGS,
               NO2_DGS, B1, 700,[])

exchange_rate = (B1.V+DGS.V+FishTank.V)*0.25/(60*60*24)

Sump_1= Sump( CO2aq_Sump, HCO3_Sump, CO32_Sump, H_Sump, OH_Sump, NH4_Sump, NH3_Sump, NO2_Sump,V_Sump, exchange_rate)
DGS = degasser( CO2aq_DGS, HCO3_DGS, CO32_DGS, H_DGS, OH_DGS, NH4_DGS, NH3_DGS, NO2_DGS, B1, 700,Sump_1)
dosing_pump_OH=Dosing_pump(OH_dosing,B1,1000)
dosing_pump_HCO3=Dosing_pump_HCO3(HCO3_dosing,B1,1000)
B1.dosing_OH=dosing_pump_OH
B1.dosing_HCO3=dosing_pump_HCO3
FishTank.Tminus = DGS

S0 = [ CO2aq_FT, HCO3_FT, CO32_FT, H_FT, OH_FT, NH4_FT, NH3_FT,NO2_FT, Fishweight_FT,Fish_Biomass_FT, CO2aq_B1, HCO3_B1, CO32_B1, H_B1, OH_B1, NH4_B1, NH3_B1, NO2_B1, AOB_B1, NOB_B1, CO2aq_DGS, HCO3_DGS, CO32_DGS, H_DGS, OH_DGS,
     NH4_DGS, NH3_DGS,NO2_DGS]


for cycle in range(0,10): #  solve the system for 10 cycles of 14 days, after each cycle the fish biomass is reduced back to 50kg/m3, the simulation is in seconds

    tspan = [14 * 24 * 60 * 60 * cycle, 14 * 24 * 60 * 60 * (cycle + 1)]

    MAlocal = solve_ivp(fun=chem, t_span=tspan, y0=S0,
                        args=(params, FishTank, B1, DGS), method="BDF")
    S0 = MAlocal.y[:, -1].tolist() #reintialize S0 with the last values of the previous cycle

    if cycle == 0:
        MA = MAlocal
    else:
        MA.y = np.append(MA.y, MAlocal.y, axis=1)
        MA.t = np.append(MA.t, MAlocal.t)

    n_out = Biomass(FishTank, MAlocal.t[-1])
    FishTank.Fish_number = FishTank.Fish_number - n_out # update the fish number
    S0[9] = FishTank.Fish_number * FishTank.Fish_weight # update the fish biomass

#plot the results
plt.plot(MA.t,(MA.y[11]*50.04+MA.y[12]*50.04+MA.y[14]*50.04),label='Alkalinity B1')
plt.plot(MA.t,(MA.y[1]*50.04+MA.y[2]*50.04+MA.y[4]*50.04),label='Alkalinity Tank')
plt.plot(MA.t,(MA.y[21]*50.04+MA.y[22]*50.04+MA.y[24]*50.04),label='Alkalinity DGS')
plt.legend()
plt.show()

# results can be store in a dataframe for further and easier exploration (e.g. TAN,TIC, CO2_removal... )
df_results = pd.DataFrame(MA.y.T, columns=['CO2aq_FT', 'HCO3_FT', 'CO32_FT', 'H_FT', 'OH_FT', 'NH4_FT', 'NH3_FT',
                                           'NO2_FT', 'Fishweight_FT', 'Fish_Biomass_FT', 'CO2aq_B1', 'HCO3_B1',
                                           'CO32_B1', 'H_B1', 'OH_B1', 'NH4_B1', 'NH3_B1', 'NO2_B1', 'AOB_B1',
                                           'NOB_B1', 'CO2aq_DGS', 'HCO3_DGS', 'CO32_DGS', 'H_DGS', 'OH_DGS',
                                           'NH4_DGS', 'NH3_DGS', 'NO2_DGS'])

# Add the time data as the first column

df_results.insert(0, 'Time', MA.t)
df_results.plot(x='Time', y=[ 'NH4_FT', 'NH4_DGS'])
plt.legend()
plt.show()
df_results.plot(x='Time', y=[ 'NO2_FT', 'NO2_DGS'])
plt.legend()
plt.show()
df_results.plot(x='Time', y=[ 'H_FT', 'H_B1', 'H_DGS'])
plt.legend()
plt.show()
plt.plot(df_results['Time'],-np.log10(df_results['H_B1']*10**-3))
plt.show()
#csv_file_path = ''
#df_results.to_csv(csv_file_path, index=False)

