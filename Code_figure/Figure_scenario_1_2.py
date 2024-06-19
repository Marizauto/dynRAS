import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

naoh_70_df = pd.read_csv('Simulation_result/simulation_results_NaOH_70.csv')
naoh_200_df= pd.read_csv('Simulation_result/simulation_results_NaOH_200.csv')
hco3_70_df = pd.read_csv('Simulation_result/simulation_results_HCO3_70.csv')
hco3_200_df  = pd.read_csv('Simulation_result/simulation_results_HCO3_200.csv')


def find_daily_min_max(df):
    df['Day'] = (df['Time']/(60*60*24)).astype(int)  
    daily_min = df.groupby('Day')['CO2aq_FT'].idxmin() 
    daily_max = df.groupby('Day')['CO2aq_FT'].idxmax()  
    return df.loc[daily_min], df.loc[daily_max]



hco3_70_daily_min, hco3_70_daily_max = find_daily_min_max(hco3_70_df)
hco3_200_daily_min, hco3_200_daily_max = find_daily_min_max(hco3_200_df)
naoh70_daily_min, naoh70_daily_max = find_daily_min_max(naoh_70_df )
naoh_200_daily_min, naoh_200_daily_max = find_daily_min_max(naoh_200_df)

height_in =10 
width_in = 8 
fig, axs = plt.subplots(4,1,figsize=(width_in,height_in))
tick_positions = np.arange(0, 25, 5)
tick_positions2= np.arange(6, 9, 0.5)

axs[0].plot(hco3_70_df['Time']/(60*60*24), hco3_70_df['CO2aq_FT']*44, color='blue', label='$CO_2$ alkalinity 70', alpha=0.7)
axs[0].plot(hco3_200_df['Time']/(60*60*24), hco3_200_df['CO2aq_FT']*44, color='black', label='$CO_2$ alkalinity 200', alpha=0.7)
axs[0].scatter(hco3_70_daily_min['Time']/(60*60*24), hco3_70_daily_min['CO2aq_FT']*44, color='blue', marker='o', s=4.5, label='Daily $CO_2$ Min/max, alkalinity 70')
axs[0].scatter(hco3_70_daily_max['Time']/(60*60*24), hco3_70_daily_max['CO2aq_FT']*44, color='blue', marker='o', s=4.5)
axs[0].scatter(hco3_200_daily_min['Time']/(60*60*24), hco3_200_daily_min['CO2aq_FT']*44, color='black', marker='o', s=4.5, label='Daily $CO_2$ Min/max, alkalinity 200')
axs[0].scatter(hco3_200_daily_max['Time']/(60*60*24), hco3_200_daily_max['CO2aq_FT']*44, color='black', marker='o', s=4.5)

axs[0].set_ylabel(r'$\mathbf{CO_2\ (mg\cdot l^{-1})}$', fontsize=12, fontweight='bold')
axs[0].set_title('a.', loc='left', fontsize=12, fontweight='bold')
axs[0].legend()
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].legend(loc='lower center', ncol=2, columnspacing=0.5, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.9),fontsize=11)

axs[0].set_ylim(0, 25)
axs[0].set_yticks(tick_positions)

axs[1].plot(naoh_70_df['Time']/(60*60*24), naoh_70_df['CO2aq_FT']*44, color='blue', label='$CO_2$ alkalinity 70', alpha=0.7)
axs[1].plot(naoh_200_df['Time']/(60*60*24), naoh_200_df['CO2aq_FT']*44, color='black', label='$CO_2$ alkalinity 200', alpha=0.7)
axs[1].scatter(naoh70_daily_min['Time']/(60*60*24), naoh70_daily_min['CO2aq_FT']*44, color='blue', marker='o', s=4.5, label='Daily $CO_2$ Min/max alkalinity 70')
axs[1].scatter(naoh70_daily_max['Time']/(60*60*24), naoh70_daily_max['CO2aq_FT']*44, color='blue', marker='o', s=4.5)
axs[1].scatter(naoh_200_daily_min['Time']/(60*60*24), naoh_200_daily_min['CO2aq_FT']*44, color='black', marker='o', s=4.5, label='Daily $CO_2$ Min/max alkalinity 200')
axs[1].scatter(naoh_200_daily_max['Time']/(60*60*24), naoh_200_daily_max['CO2aq_FT']*44, color='black', marker='o', s=4.5)

axs[1].set_ylabel(r'$\mathbf{CO_2\ (mg\cdot l^{-1})}$', fontsize=12, fontweight='bold')
axs[1].set_ylim(0, 25)
axs[1].set_xlim(0, 141)
axs[1].set_yticks(tick_positions)
axs[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs[1].legend(loc='lower center', ncol=2,columnspacing=0.5, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.9),fontsize=11)


axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)
axs[2].set_title('c.', loc='left', fontsize=12, fontweight='bold')
axs[2].plot(hco3_70_df['Time']/(60*60*24), -np.log10(hco3_70_df['H_FT']*10**-3), color='blue', label='pH alkalinity 70', alpha=0.8)
axs[2].plot(hco3_200_df['Time']/(60*60*24), -np.log10(hco3_200_df['H_FT']*10**-3), color='black', label='pH alkalinity 200', alpha=0.8)
axs[2].set_ylim(6, 8)
axs[3].set_ylabel('pH', fontsize=12, fontweight='bold')
axs[3].set_title('d.', loc='left', fontsize=12, fontweight='bold')
axs[3].plot(naoh_70_df['Time']/(60*60*24), -np.log10(naoh_70_df['H_FT']*10**-3), color='blue', label='pH alkalinity 70', alpha=0.8)
axs[3].plot(naoh_200_df['Time']/(60*60*24), -np.log10(naoh_200_df['H_FT']*10**-3), color='black', label='pH alkalinity 200', alpha=0.8)
axs[2].set_ylabel('pH', fontsize=12, fontweight='bold')
axs[3].spines['top'].set_visible(False)
axs[3].spines['right'].set_visible(False)
axs[2].spines['top'].set_visible(False)
axs[2].spines['right'].set_visible(False)
axs[3].set_ylim(6, 8)
axs[3].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs[2].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94),fontsize=11)
axs[3].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94),fontsize=11)
axs[0].set_xlim(0, 141)
axs[2].set_xlim(0, 141)
axs[3].set_xlim(0, 141)
axs[2].set_yticks(tick_positions2)
axs[3].set_yticks(tick_positions2)
plt.tight_layout()
fig.tight_layout(pad=1.0)
plt.savefig('figure/scenario_1_2/Comparison_HCO3_NaOH.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()

import matplotlib.lines as mlines


height_in =10 
width_in = 8  
fig, axs = plt.subplots(4,1,figsize=(width_in,height_in))
axs = axs.flatten()
labels = ['a.', 'b.', 'c.', 'd.']
titles = ['$HCO_3^-$ (alkalinity 70 mg.$L^{-1}$)', '$HCO_3^-$ (alkalinity 200 mg.$L^{-1}$)', 'NaOH (alkalinity 70 mg.$L^{-1}$)', 'NaOH (alkalinity 200 mg.$L^{-1}$)']
datasets = [naoh_70_df, naoh_200_df, hco3_70_df, hco3_200_df]
time_in_days = lambda time: time / (60 * 60 * 24)
for ax, df, title, label in zip(axs, datasets, titles, labels):
    ax.plot(time_in_days(df['Time']), (df['NH4_FT'] + df['NH3_FT'])*14, color='black', linewidth=1.5)
    ax.set_ylabel('TAN $\mathbf{(mg\cdot L^{-1})}$', color='black', fontweight='bold')
    ax.tick_params(axis='y', labelcolor='black')
    ax.set_title(title)
    ax.set_ylim(0, 3)
    ax.set_xlim(0,141)
    ax2 = ax.twinx()
    ax2.plot(time_in_days(df['Time']), df['NH3_FT']*14*1e3, color='blue', linestyle='-', linewidth=1, alpha=0.7)
    ax2.set_ylabel('NH3-N $\mathbf{(μg\cdot L^{-1})}$', color='blue', fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='blue')
    ax2.set_ylim(0, 25)
    ax.text(-0.1, 1.1, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax.set_xlabel('Time (days)', fontweight='bold')
black_line = mlines.Line2D([], [], color='black', linewidth=1.5, label='TAN ($mg\cdot L^{-1}$)')
blue_line = mlines.Line2D([], [], color='blue', linewidth=1, linestyle='-', alpha=0.7, label='NH3-N $\mathbf{(μg\cdot L^{-1})}$')
fig.legend(handles=[black_line, blue_line], loc='lower center', ncol=2, bbox_to_anchor=(0.5, 0))

fig.tight_layout()
plt.subplots_adjust(top=0.9, bottom=0.15)
plt.savefig('figure/scenario_1_2/TAN_comparison_OH_HCO3.png', dpi=300, bbox_inches='tight',pad_inches=0)

naoh_70_df['TIC_removal'] = (((naoh_70_df['CO2aq_FT']+naoh_70_df['HCO3_FT']+naoh_70_df['CO32_FT'])-(naoh_70_df['CO2aq_DGS']+naoh_70_df['HCO3_DGS']+naoh_70_df['CO32_DGS']))/(naoh_70_df['CO2aq_FT']+naoh_70_df['HCO3_FT']+naoh_70_df['CO32_FT']))*100
naoh_200_df['TIC_removal'] = (((naoh_200_df['CO2aq_FT']+naoh_200_df['HCO3_FT']+naoh_200_df['CO32_FT'])-(naoh_200_df['CO2aq_DGS']+naoh_200_df['HCO3_DGS']+naoh_200_df['CO32_DGS']))/(naoh_200_df['CO2aq_FT']+naoh_200_df['HCO3_FT']+naoh_200_df['CO32_FT']))*100
naoh_70_df['CO2_removal'] = (((naoh_70_df['CO2aq_FT'])-(naoh_70_df['CO2aq_DGS']))/(naoh_70_df['CO2aq_FT']))*100
naoh_200_df['CO2_removal'] = (((naoh_200_df['CO2aq_FT'])-(naoh_200_df['CO2aq_DGS']))/(naoh_200_df['CO2aq_FT']))*100
hco3_70_df['TIC_removal'] = (((hco3_70_df['CO2aq_FT']+hco3_70_df['HCO3_FT']+hco3_70_df['CO32_FT'])-(hco3_70_df['CO2aq_DGS']+hco3_70_df['HCO3_DGS']+hco3_70_df['CO32_DGS']))/(hco3_70_df['CO2aq_FT']+hco3_70_df['HCO3_FT']+hco3_70_df['CO32_FT']))*100
hco3_200_df['TIC_removal'] = (((hco3_200_df['CO2aq_FT']+hco3_200_df['HCO3_FT']+hco3_200_df['CO32_FT'])-(hco3_200_df['CO2aq_DGS']+hco3_200_df['HCO3_DGS']+hco3_200_df['CO32_DGS']))/(hco3_200_df['CO2aq_FT']+hco3_200_df['HCO3_FT']+hco3_200_df['CO32_FT']))*100
hco3_70_df['CO2_removal'] = (((hco3_70_df['CO2aq_FT'])-(hco3_70_df['CO2aq_DGS']))/(hco3_70_df['CO2aq_FT']))*100
hco3_200_df['CO2_removal'] = (((hco3_200_df['CO2aq_FT'])-(hco3_200_df['CO2aq_DGS']))/(hco3_200_df['CO2aq_FT']))*100
def calculate_stats(df):
    return pd.Series({
        'Mean_TIC_removal': df['TIC_removal'].mean(),
        'Std_TIC_removal': df['TIC_removal'].std(),
        'Mean_CO2_removal': df['CO2_removal'].mean(),
        'Std_CO2_removal': df['CO2_removal'].std()
    })
stats_HCO3_70 = calculate_stats(hco3_70_df)
stats_HCO3_200 = calculate_stats(hco3_200_df)
summary_df_HCO3 = pd.DataFrame({
    'Concentration': [70, 200],
    'Mean_TIC_removal': [stats_HCO3_70['Mean_TIC_removal'], stats_HCO3_200['Mean_TIC_removal']],
    'Std_TIC_removal': [
                        stats_HCO3_70['Std_TIC_removal'], stats_HCO3_200['Std_TIC_removal']],
    'Mean_CO2_removal': [
                         stats_HCO3_70['Mean_CO2_removal'],  stats_HCO3_200['Mean_CO2_removal']],
    'Std_CO2_removal': [
                        stats_HCO3_70['Std_CO2_removal'], stats_HCO3_200['Std_CO2_removal']]
})
stats_NaOH_70 = calculate_stats(naoh_70_df)
stats_NaOH_200 = calculate_stats(naoh_200_df)
summary_df_NaOH = pd.DataFrame({
    'Concentration': [70, 200],
    'Mean_TIC_removal': [stats_NaOH_70['Mean_TIC_removal'], stats_NaOH_200['Mean_TIC_removal']],
    'Std_TIC_removal': [
                        stats_NaOH_70['Std_TIC_removal'],  stats_NaOH_200['Std_TIC_removal']],
    'Mean_CO2_removal': [
                         stats_NaOH_70['Mean_CO2_removal'], stats_NaOH_200['Mean_CO2_removal']],
    'Std_CO2_removal': [
                        stats_NaOH_70['Std_CO2_removal'], stats_NaOH_200['Std_CO2_removal']]
})

width = 0.35

fig, ax = plt.subplots(figsize=(10, 7))

positions_HCO3 = np.arange(len(summary_df_HCO3)) - width/2


positions_NaOH = np.arange(len(summary_df_NaOH)) + width/2


co2_bars_HCO3 = ax.bar(positions_HCO3 + width/4, summary_df_HCO3['Mean_CO2_removal'], width/2, yerr=summary_df_HCO3['Std_CO2_removal'], capsize=5, label='CO2 Removal - HCO3', color='darkcyan', ecolor='black',edgecolor='black')
co2_bars_NaOH = ax.bar(positions_NaOH + width/4, summary_df_NaOH['Mean_CO2_removal'], width/2, yerr=summary_df_NaOH['Std_CO2_removal'], capsize=5, label='CO2 Removal - NaOH', color='gray', ecolor='black',edgecolor='black')

tic_bars_NaOH = ax.bar(positions_NaOH - width/4, summary_df_NaOH['Mean_TIC_removal'], width/2, yerr=summary_df_NaOH['Std_TIC_removal'], capsize=5, label='TIC Removal - NaOH', color='gray', hatch='/',ecolor='black',edgecolor='black')
tic_bars_HCO3 = ax.bar(positions_HCO3 - width/4, summary_df_HCO3['Mean_TIC_removal'], width/2, yerr=summary_df_HCO3['Std_TIC_removal'], capsize=5, label='TIC Removal - HCO3', color='darkcyan', hatch='\\', ecolor='black',edgecolor='black')


ax.set_ylabel('Removal Efficiency (%)', fontweight='bold')
ax.set_xticks(np.arange(len(summary_df_HCO3)))
ax.set_xticklabels(summary_df_HCO3['Concentration'])
plt.xlabel(r'$\mathbf{Alkalinity\ (mg\cdot L^{-1})\ as\ CaCO_3 }$', fontweight='bold')
ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.25))
ax.set_axisbelow(True)
ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('figure/scenario_1_2/HCO3_NaOH_TIC.png', dpi=300)
plt.show()
plt.show()
