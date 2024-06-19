# ![image](https://github.com/Marizauto/dynRAS/assets/128140640/5b76d9d7-d748-4974-bbe4-4270cd8e15a6)



The dynRAS model is a model developed to simulate recirculating aquaculture systems (RAS) with a particular emphasis on alkalinity and pH control. It is build on prior research effort (Pedersen 2018 and Wik 2008) with new functionalities, including the incorporation of the carbonate system to accurately represent alkalinity and pH fluctuations. Developed within the scope of the RASTOOLS project, the model aims to investigate the dynamic of RAS. The RAS components included in dynRAS are simplified to the fish tank, biofilter, and degasser, based on an experimental set-up by Jafari et al. (2024). The model operates on Python 3.8 and necessitates the scipy, numpy, and math libraries.

* Python scripts entitled "Class_definition_and_solver.py" is the main script to run to solve the system of differential equations describing the change of the system over time. 
* ChemODE_BIO.py,ChemODE_DGS.py,ChemODE_Fish.py contains the system of differential equation for the biofilter, degasser and fish tank. ChemODE_BIO.py can be modify to simulate different dosing scenario.
* Fish_growth.py and Growth_Bacteria.py include the EDO for the fish growth and the growth of the bacteria. Only AOB and NOB are included in the model.
* TAN_prod.py contains a function to simulate TAN excretion from fish
* Biomass_function.py simulates fish removal using a threshold for the maximum biomass present in the system. The function is called every 14 days.
* params.py is the list of paramters for the model

While fully functional, the model still needs improvement. Notably, numerical issues can be encounter for particular scenarios, particularly when controling the pH through HCO3 for this scenario, the reduction of the steepness for the answer of the dosage system seems to fix these numerical issue however the dosage system becomes less efficient. The chemical equilibriums have to be consider when setting the innitial conditions, if not the solver will struggle to converge on a solution. e.g. Setting a very high pH, a very high $CO_2$ and a very low $HCO_3^-$ will lead to numerical instabilities.

How to simulate scenarios from the paper ?

* Simulation of the experimental set-up from Jaffari et al. (2024)
When running the "Class_definition_and_solver.py" file the scenario by default is the experimental set-up from Jaffari et. al.(2024). The 
* Simulating scenario 1
To run scenario 1, in "Class_definition_and_solver.py" comment line 94 "dY = chemODE_BIO(self, params,t)" and uncomment line dY = chemODE_BIO_HCO3(self, params, t). This will run scenario 1 for alkalinity 200. To run scenario 1 at alkalinity 70, in the file ChemODE_BIO.py uncomment line 137 and comment line 136. Then go to the python file "Class_definition_and_solver.py", comment line 168,187,200 and uncomment line 169,188,201. Run the python file "Class_definition_and_solver.py" once again the out
