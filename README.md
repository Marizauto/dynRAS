# ![image](https://github.com/Marizauto/dynRAS/assets/128140640/5b76d9d7-d748-4974-bbe4-4270cd8e15a6)

The dynRAS model is a model developed to simulate recirculating aquaculture systems (RAS) with a particular emphasis on alkalinity and pH control. It is build on prior research effort (Pedersen 2018 and Wik 2008) with new functionalities, including the incorporation of the carbonate system to accurately represent alkalinity and pH fluctuations. Developed within the scope of the RASTOOLS project, the model aims to investigate the dynamic of RAS. The RAS components included in dynRAS are simplified to the fish tank, biofilter, and degasser, based on an experimental set-up by Jafari et al. (2024). The model operates on Python 3.8 and necessitates the scipy, numpy, pandas, math and matplotlib libraries.The model operates on Python 3.8 and requires the `scipy`, `numpy`, `pandas`, `math`, and `matplotlib` libraries.

## Model Components and Files

The model can be found within the folder `DynRAS_model` and includes the following files:

- **Class_definition_and_solver.py**: The main script to solve the system of differential equations describing the changes in the system over time.
- **ChemODE_BIO.py, ChemODE_DGS.py, ChemODE_Fish.py**: Contain the system of differential equations for the biofilter, degasser, and fish tank, respectively. `ChemODE_BIO.py` can be modified to simulate different dosing scenarios.
- **Fish_growth.py and Growth_Bacteria.py**: Include the differential equations (ODE) for the growth of fish and bacteria (AOB and NOB).
- **TAN_prod.py**: Contains a function to simulate Total Ammonia Nitrogen (TAN) excretion from fish.
- **Biomass_function.py**: Simulates fish removal using a threshold for the maximum biomass present in the system. This function is called every 14 days.
- **params.py**: Lists the parameters for the model.

## Running the Model

### Default Simulation

To simulate the experimental set-up from Jafari et al. (2024), simply run the `Class_definition_and_solver.py` file. The default scenario is set to match the experimental conditions.

### Scenario 1

1. **For Alkalinity 200:**
   - In `Class_definition_and_solver.py`, comment out lines 94, 96, 97, 98; uncomment line 95, `dY = chemODE_BIO_HCO3(self, params, t)`, uncomment out lines 168, 187, 200, comment out lines 169, 188, 201.
   - In `ChemODE_BIO.py`, comment out line 107 and uncomment line 106.
   - Run `Class_definition_and_solver.py`.

2. **For Alkalinity 70:**
   - In `ChemODE_BIO.py`, uncomment line 107 and comment line 106.
   - In `Class_definition_and_solver.py`, comment out lines 168, 187, 200, uncomment lines 169, 188, 201.
   - Run `Class_definition_and_solver.py`.

### Scenario 2

1. **For Alkalinity 200:**
   - In `Class_definition_and_solver.py`, comment out lines 94, 95, 97, 98,169, 188, 201; uncomment line 96 `dY = chemODE_BIO_NaOH(self, params, t)` and lines 168, 187, 200.
   - In `ChemODE_BIO.py`, uncomment line 136 and comment line 137.
   - Run `Class_definition_and_solver.py`.

3. **For Alkalinity 70:**
   - In `ChemODE_BIO.py`, uncomment line 137 and comment line 136.
   - In `Class_definition_and_solver.py`, comment out lines 168, 187, 200; uncomment lines 169, 188, 201.
   - Run `Class_definition_and_solver.py`.

### Scenario 3

1. **For Alkalinity 200:**
   - In `Class_definition_and_solver.py`, comment out lines 94, 95, 96, 98,169, 188, 201; uncomment line 97 `dY = chemODE_BIO_alk_control(self, params, t)` and lines 168, 187, 200.
   - In `ChemODE_BIO.py`, uncomment line 168 and comment line 169.
   - Run `Class_definition_and_solver.py`.

2. **For Alkalinity 70:**
   - In `ChemODE_BIO.py`, comment line 168 and uncomment line 169.
   - In `Class_definition_and_solver.py`, comment out lines 168, 187, 200; uncomment lines 169, 188, 201.
   - Run `Class_definition_and_solver.py`.

### pH Control Instead of Alkalinity Control

   - In `Class_definition_and_solver.py`, uncomment line 98 and comment lines 94, 95, 96, 97 to switch to pH control.
   - In `ChemODE_BIO.py`, to control pH with $HCO_3^-$, uncomment lines 204 and 211, and comment lines 212 and 206.
   - To control pH with NaOH, uncomment lines 212 and 206, and comment lines 204 and 211.
   - Adjust the pH threshold on line 202 if necessary.
   - Initial conditions in `Class_definition_and_solver.py` might need adaptation.

When running the Class_definition_and_solver.py script, it generates two datafram that will be saved in two csv file 'results.csv' and 'dosing.csv'. The 'results.csv' file contains the simulation results for each species in each compartment of the system, first column is the time ('Time') and the name of the other columns can be found from line 267 to line 271. The name of the columns is build with the name of the species first, followed by the compartment '_FT' for the fish tank, '_B1' for the biofilter and '_DGS' for the Degasser. (e.g. 'CO2aq_FT' is the $CO_2$ result in the fish tank). Time is in second, the chemical species concentration are in mmol.L^{-1} you can multiply the species concentration by their respective molecular weight given from line 274 to line 279 to get the concentration in mg.l^{-1}. The 'dosing.csv' contains 3 columns with the first column containing  the simulation time, second column the dosing rate of OH and the dosing rate of HCO3 in mmol.l^{-1}.s^{-1}. To prevent overwriting this file with each run, change the filenames on line 281 and 314 of the script. This will ensure that each simulation's results are saved to a new file. In addition, from line 288 you have exemple on how to plot the results from the model. The simulation will automatically plot results for $CO_2$, alkalinity and TAN in mg per litre from the fish tank as well as the dosing rate of the two supplements over the simulation time.
## Generating Figures

### Data

To generate the same figures as in the paper, use either the output from the model simulations or the data available in the `Code_figure` directory. This directory contains:
- `Experimental_data`: Folder with experimental data collected.
- `Simulation_result`: Folder with simulation outputs for each scenario presented in the paper.

### Creating Figures

The scripts used to generate the figures are in the `Code_figure` directory:
- `Jafari_et_al_2024.py`: Generates figures for model validation.
- `Figure_scenario_1_2.py`: Generates figures for scenarios 1 and 2.
- `Scenario_3.py`: Generates figures for scenario 3.

Running each file will create the figures and store them in a directory named `figure`, which contains subfolders for the figures corresponding to each script.

While fully functional, the model still needs improvement. Notably, numerical issues can be encounter for particular scenarios, particularly when controling the pH through HCO3 for this scenario, the reduction of the steepness for the answer of the dosage system seems to fix these numerical issue however the dosage system becomes less efficient. The chemical equilibriums have to be consider when setting the innitial conditions, if not the solver will struggle to converge on a solution. e.g. Setting a very high pH, a very high $CO_2$ and a very low  $HCO_3^-$ will lead to numerical instabilities.
