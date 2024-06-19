import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12

alk_70_df_dosage = pd.read_csv('Simulation_result/dosing_pH_Alk_control_70.csv')
alk_200_df_dosage= pd.read_csv('Simulation_result/dosing_pH_Alk_control_200.csv')
alk_70_df=pd.read_csv('Simulation_result/simulation_results_pH_Alk_control_70.csv')
alk_200_df=pd.read_csv('Simulation_result/simulation_results_pH_Alk_control_200.csv')
print(alk_70_df_dosage.head())
def find_daily_min_max(df):
    df['Day'] = (df['Time']/(60*60*24)).astype(int)
    daily_min = df.groupby('Day')['CO2aq_FT'].idxmin()
    daily_max = df.groupby('Day')['CO2aq_FT'].idxmax()
    return df.loc[daily_min], df.loc[daily_max]



df_70_daily_min, df_70_daily_max = find_daily_min_max(alk_70_df)
df_200_daily_min, df_200_daily_max = find_daily_min_max(alk_200_df)

height_in =10/2
width_in = 8
fig, axs = plt.subplots(2, 1, figsize=(width_in,height_in))


axs[0].plot(alk_200_df_dosage['Time']/(60*60*24), alk_200_df_dosage['OH_dosing']*10**3, color='black', label='OH dosing')
axs[0].plot(alk_200_df_dosage['Time']/(60*60*24), alk_200_df_dosage['HCO3_dosing']*10**3, color='blue', label='HCO3 dosing')
axs[0].set_ylabel(r'$\mathbf{Dosing \ rate}$' + '\n' + r'$\mathbf{(\mu mol\cdot L^{-1}\cdot s^{-1})}$', fontsize=12, fontweight='bold')
axs[0].set_title('a.', loc='left', fontsize=12, fontweight='bold')
axs[0].set_ylim(0,0.22)
axs[0].set_title('Alkalinity 200', loc='center', fontsize=12, fontweight='bold')
axs[0].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.8),fontsize=11)

axs[1].plot(alk_70_df_dosage['Time']/(60*60*24), alk_70_df_dosage['OH_dosing']*10**3, color='black', label='OH dosing')
axs[1].plot(alk_70_df_dosage['Time']/(60*60*24), alk_70_df_dosage['HCO3_dosing']*10**3, color='blue', label='HCO3 dosing')
axs[1].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs[1].set_ylabel(r'$\mathbf{Dosing \ rate}$' + '\n' + r'$\mathbf{(\mu mol\cdot L^{-1}\cdot s^{-1})}$', fontsize=12, fontweight='bold')
axs[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs[1].set_title('Alkalinity 70', loc='center', fontsize=12, fontweight='bold')
axs[1].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.8),fontsize=11)
axs[1].set_ylim(0,0.22)
axs[0].set_xlim(0, 141)
axs[1].set_xlim(0, 141)
fig.tight_layout(pad=1.0)
plt.savefig('figure/scenario_3/dosing.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()

fig, axs = plt.subplots(2,1, figsize=(width_in,height_in))
axs[0].scatter(df_70_daily_min['Time']/(60*60*24), df_70_daily_min['CO2aq_FT']*44, color='blue', marker='o', s=4.5, label='Daily $CO_2$ Min/max, alkalinity 70')
axs[0,].scatter(df_70_daily_max['Time']/(60*60*24), df_70_daily_max['CO2aq_FT']*44, color='blue', marker='o', s=4.5)
axs[0].scatter(df_200_daily_min['Time']/(60*60*24), df_200_daily_min['CO2aq_FT']*44, color='black', marker='o', s=4.5, label='Daily $CO_2$ Min/max, alkalinity 200')
axs[0].scatter(df_200_daily_max['Time']/(60*60*24), df_200_daily_max['CO2aq_FT']*44, color='black', marker='o', s=4.5)
axs[0].plot(alk_70_df['Time']/(60*60*24), alk_70_df['CO2aq_FT']*44, color='blue', label='alkalinity 70',linewidth=1.5)
axs[0].plot(alk_200_df['Time']/(60*60*24), alk_200_df['CO2aq_FT']*44, color='black', label='alkalinity 200',alpha=0.8,linewidth=1)
axs[0].set_ylabel(r'$\mathbf{CO_2 \ (mg\cdot L^{-1})}$', fontsize=12, fontweight='bold')
axs[0].set_title('a.', loc='left', fontsize=12, fontweight='bold')
axs[0].legend()
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].set_ylim(0,25)
axs[0].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.6),fontsize=11)

axs[1].plot(alk_70_df['Time']/(60*60*24), -np.log10(alk_70_df['H_FT']*10**-3), color='blue', label='alkalinity 70')
axs[1].plot(alk_200_df['Time']/(60*60*24), -np.log10(alk_200_df['H_FT']*10**-3), color='black', label='alkalinity 200')
axs[1].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs[1].set_ylabel(r'$\mathbf{pH}$', fontsize=12, fontweight='bold')
axs[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs[1].set_ylim(6,8)
axs[1].legend()
axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)
axs[1].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.8),fontsize=11)
axs[0].set_xlim(0, 141)
axs[1].set_xlim(0, 141)
fig.tight_layout(pad=1.0)
plt.savefig('figure/scenario_3/CO2_pH.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()

fig, axs = plt.subplots(2,1, figsize=(width_in,height_in))
axs[0].plot(alk_200_df['Time']/(60*60*24),(alk_200_df['NH4_FT'] + alk_200_df['NH3_FT'])*12, color='black',linewidth=1.5,)
axs[0].set_ylabel(r'$\mathbf{TAN\ (mg\cdot L^{-1}\cdot s^{-1})}$', fontsize=12, fontweight='bold')
axs[0].set_title('a.', loc='left', fontsize=12, fontweight='bold')
axs[0].set_ylim(0,3)
axs[0].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94),fontsize=11)
ax2 = axs[0].twinx()
ax2.plot(alk_200_df['Time']/(60*60*24), alk_200_df['NH3_FT']*12*1e3, color='blue',linewidth=1,alpha=0.7)
ax2.set_ylabel(r'$\mathbf{NH_3-N \ (\mu g\cdot L^{-1})}$', fontsize=12, fontweight='bold',color='blue')
ax2.set_ylim(0,25)
axs[0].set_title('Alkalinity 200', loc='center', fontsize=12, fontweight='bold')
axs[1].plot(alk_70_df['Time']/(60*60*24),(alk_70_df['NH4_FT'] + alk_70_df['NH3_FT'])*12, color='black',linewidth=1.5)
axs[1].set_xlabel('Time (days)', fontsize=12, fontweight='bold')
axs[1].set_ylabel(r'$\mathbf{TAN\ (mg\cdot L^{-1}\cdot s^{-1})}$', fontsize=12, fontweight='bold')
axs[1].set_title('b.', loc='left', fontsize=12, fontweight='bold')
axs[1].set_ylim(0,3)
axs[1].legend(loc='lower center', ncol=2, framealpha=0, facecolor='white', bbox_to_anchor=(0.5, 0.94),fontsize=11)
ax2 = axs[1].twinx()
ax2.plot(alk_70_df['Time']/(60*60*24), alk_70_df['NH3_FT']*12*1e3,  color='blue',linewidth=1,alpha=0.7)
ax2.set_ylabel(r'$\mathbf{NH_3-N \ (\mu g\cdot L^{-1})}$', fontsize=12, fontweight='bold',color='blue')
ax2.set_ylim(0,25)
axs[1].set_title('Alkalinity 70', loc='center', fontsize=12, fontweight='bold')
plt.tight_layout()
axs[0].set_xlim(0, 141)
axs[1].set_xlim(0, 141)
plt.savefig('figure/scenario_3/TAN_NH3.png', dpi=300, bbox_inches='tight',pad_inches=0)
plt.show()

