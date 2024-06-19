import numpy as np

def Biomass(FishTank,Y):


    n_out=0

    #if round(t)%1209600==0:

    if FishTank.Fish_Biomass> FishTank.Fish_Biomass_Max:
        dt=FishTank.Fish_Biomass-FishTank.Fish_Biomass_Max
        n_out = round(dt/FishTank.Fish_weight)
        #FishTank.Fish_number = FishTank.Fish_number-n_out
        #Fish_Biomass=FishTank.Fish_number*FishTank.Fish_weight

    return n_out