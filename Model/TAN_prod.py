import math
from scipy.stats import beta


def NH4_rate(self, t):
    day = t / (60 * 60 * 24)  # Convert seconds to days
    interval = 24 * (day - math.floor(day))  # Current hour in the day

    # Total TAN excretion
    Feed_total = (0.027* self.Fish_Biomass) #g

    TAN_total = (((Feed_total * 0.46 * 0.092) / 14)) #mol
    TAN_total= TAN_total*1000 #convert to mmol
    # Beta distribution parameters
    alpha_param, beta_param = 3,3 # Adjust these values as needed
    proportion = beta.pdf(((interval)/ 24), alpha_param, beta_param)
    min_proportion = 0.1 #ensure a non zero value for TAN excretion
    modified_proportion = proportion * (1 - min_proportion) + min_proportion

    return modified_proportion * TAN_total