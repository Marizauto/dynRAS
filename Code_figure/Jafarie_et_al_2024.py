
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

data = pd.read_csv('Experimental_data/co2_alkalinity.csv')
simulation = pd.read_csv('Simulation_result/Jafarie_et_al_2024.csv')
data2 = pd.read_csv('Experimental_data/ph_data.csv')

height_in =10-(10/3)
width_in = 8
yticks = np.arange(0, 25, 5)
yticks2 = np.arange(6, 9, 0.5)
yticks3 = np.arange(0, 500, 100)
data_weight = {
    'bulk_weight': ["27/09/22", "12/10/22", "28/10/22", "11/11/22", "26/11/22", "09/12/22", "20/12/22", "05/01/23","17/01/23"],
    'nb_fish': [352, 241, 198, 154, 126, 110, 93, 83, 0]
}
df = pd.DataFrame(data_weight)
df['bulk_weight'] = pd.to_datetime(df['bulk_weight'], format='%d/%m/%y')

experiment_start = pd.to_datetime("29/08/22 08:00:00", format='%d/%m/%y %H:%M:%S')

df['days_since_start'] = (df['bulk_weight'] - experiment_start).dt.total_seconds() / (60 * 60 * 24)

fig1, axs1 = plt.subplots(3, 1, figsize=(width_in, height_in))
for ax in axs1:
    for i, day in enumerate(range(0, 141, 14)):
        if i == 0:
            ax.axvline(x=day, color='grey', linestyle='--', linewidth=0.5, label='Fish removal simulation')
        else:
            ax.axvline(x=day, color='grey', linestyle='--', linewidth=0.5)
for ax in axs1:
    for day in df['days_since_start']:
        ax.axvline(x=day, color='red', linestyle='--', linewidth=1, label='Fish removal data' if day == df['days_since_start'].iloc[0] else "")
# Panel a
axs1[0].plot(data['time'], data['co2_mgl'], label='$CO_2$ data', color='black')
axs1[0].plot(simulation['Time']/(60*60*24), simulation['CO2aq_FT']*44, label="$CO_2$ simulation", color='blue')
axs1[0].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94))
axs1[0].set_title('a.', loc='left', fontsize=12, fontweight='bold')
axs1[0].set_ylabel(r'$\mathbf{CO_2\ (mg\cdot l^{-1})}$', fontsize=12, fontweight='bold')
axs1[0].set_ylim(0, 25)
axs1[0].set_yticks(yticks)
# Panel b
axs1[1].plot(data2['time'], data2['ph'], label='pH data', color='black')
axs1[1].plot(simulation['Time']/(60*60*24), -np.log10(simulation['H_FT']*10**-3), label="pH simulation", color='blue')
axs1[1].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94))
axs1[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs1[1].set_ylabel('pH', fontsize=12, fontweight='bold')
axs1[1].set_ylim(6, 8.5)
axs1[1].set_yticks(yticks2)

# Panel c
axs1[2].plot(data['time'], data['alkalinity_mgl'], label='Alkalinity data', color='black')
axs1[2].plot(simulation['Time']/(60*60*24), simulation['OH_FT']*50.04 + simulation['HCO3_FT']*50.04 + 2*simulation['CO32_FT']*50.04, label="Alkalinity simulation", color='blue')
axs1[2].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94))
axs1[2].set_title('c.', loc='left', fontsize=12, fontweight='bold')
axs1[2].set_ylabel(r'$\mathbf{Alkalinity} $' + '\n' + r'$\mathbf{(mg\cdot l^{-1}CaCO_3)}$', color='black',fontweight='bold')
axs1[2].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs1[2].set_ylim(0, 500)
axs1[2].set_yticks(yticks3)


axs1[0].spines['top'].set_visible(False)
axs1[0].spines['right'].set_visible(False)
axs1[1].spines['top'].set_visible(False)
axs1[1].spines['right'].set_visible(False)
axs1[2].spines['top'].set_visible(False)
axs1[2].spines['right'].set_visible(False)
axs1[0].set_xlim(0, 141)
axs1[1].set_xlim(0, 141)
axs1[2].set_xlim(0, 141)
fig1.tight_layout(pad=1.0)
plt.savefig('figure/Jafarie_et_al_settings/Comparison_simulation_data.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()

height_in =10/2
width_in = 8
xticks = np.arange(0, 21, 1)
yticks = np.arange(0, 25, 5)
fig1, axs1 = plt.subplots(2, 1, figsize=(width_in, height_in))
axs1[0].set_xlim(0, 20)
axs1[0].plot(data['time'], data['co2_mgl'], label='$CO_2$ data', color='black')
axs1[0].plot(simulation['Time']/(60*60*24), simulation['CO2aq_FT']*44, label="$CO_2$ simulation", color='blue')
axs1[0].set_title('a.',loc='left', fontsize=12, fontweight='bold')
axs1[0].set_ylabel(r'$\mathbf{CO_2\ (mg\cdot l^{-1})}$', fontsize=12, fontweight='bold')
for start in np.arange(0, 20, 1):
    axs1[0].axvspan(start, start + 0.5, color='gray', alpha=0.3)
feeding_patch = Patch(color='gray', alpha=0.3, label='Feeding period')
axs1[0].legend(handles=[feeding_patch, *axs1[0].get_legend_handles_labels()[0]], loc='lower center', ncol=3, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.95))
axs1[1].plot(data2['time'], data2['ph'], label='pH data', color='black')
axs1[1].plot(simulation['Time']/(60*60*24), -np.log10(simulation['H_FT']*10**-3), label="pH simulation", color='blue')
axs1[1].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94))
axs1[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs1[1].set_ylabel('pH', fontsize=12, fontweight='bold')
axs1[1].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs1[1].set_ylim(6, 8.5)
for start in np.arange(0, 20, 1):
    axs1[1].axvspan(start, start + 0.5, color='gray', alpha=0.3)
axs1[1].set_xlim(0, 20)
feeding_patch = Patch(color='gray', alpha=0.3, label='Feeding period')
fig1.tight_layout(pad=1.0)
axs1[1].legend(handles=[feeding_patch, *axs1[0].get_legend_handles_labels()[0]], loc='lower center', ncol=3, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.95))
axs1[0].set_xticks(xticks)
axs1[0].set_yticks(yticks)
axs1[1].set_xticks(xticks)
axs1[1].set_yticks(np.arange(6, 9, 0.5))
plt.savefig('figure/Jafarie_et_al_settings/20days.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()
xticks = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 0, 1, 2, 3, 4, 5, 6, 7]
simulation['Time_hours'] = (simulation['Time'] / 3600) % 24
simulation['NH4_NH3_Combined'] = simulation['NH4_FT']*14 + simulation['NH3_FT']*14
hourly_means_nh4_nh3 = simulation.groupby(simulation['Time_hours'].astype(int))['NH4_NH3_Combined'].mean()
hourly_stds_nh4_nh3 = simulation.groupby(simulation['Time_hours'].astype(int))['NH4_NH3_Combined'].std()


Timeseries = pd.read_csv(r"C:\Users\marie\OneDrive\Bureau\STATSfile\Timeserie_module_1.csv")
Timeseries['date'] = pd.to_datetime(Timeseries['date'])


Timeseries['T0_as_8AM_hour'] = (Timeseries['date'].dt.hour + 16) % 24
T0_as_8AM_hourly_tan_means = Timeseries.groupby('T0_as_8AM_hour')['TAN_DF'].mean()
T0_as_8AM_hourly_tan_stds = Timeseries.groupby('T0_as_8AM_hour')['TAN_DF'].std()


fig, ax = plt.subplots(figsize=(7, 4))
bar_width = 0.35
index_simulation = np.arange(24)


ax.bar(index_simulation - bar_width/2, hourly_means_nh4_nh3.reindex(index_simulation).fillna(0), bar_width, yerr=hourly_stds_nh4_nh3.reindex(index_simulation).fillna(0), capsize=5, label='Simulation', color='black', error_kw=dict(elinewidth=1, ecolor='gray'))


index_timeseries = T0_as_8AM_hourly_tan_means.index
ax.bar(index_timeseries + bar_width/2, T0_as_8AM_hourly_tan_means, bar_width, yerr=T0_as_8AM_hourly_tan_stds, capsize=5, label='Measurement', color='grey', edgecolor='black', hatch='//', error_kw=dict(elinewidth=1, ecolor='black'))
ax.set_axisbelow(True)
ax.set_xlabel('Hour of the day', fontsize=12)
ax.set_ylabel('Mean TAN Values in $mg.L^{-1}$', fontsize=12)
ax.set_xticks(index_simulation)  # Set x-ticks to cover all 24 hours
ax.set_xticklabels([f'{i}:00' for i in xticks], rotation=45)
ax.legend()
plt.grid(axis='y', linestyle='--')
plt.tight_layout()
fig.tight_layout(pad=1.0)
plt.savefig('figure/Jafarie_et_al_settings/TAN.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()