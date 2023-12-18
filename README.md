## Kinetic Monte Carlo Simulator for ABVI Binary Alloy System ##




### Introduction ###

Perform kinetic Monte Carlo simulations for a binary A-B alloy with defects of vacancies and interstitials. The atomic diffusions are assumed via defect mechanisms including: defect jumps, recombinations, annihilations at defect sinks. Frenkel-Pair generation is implemented for defect insertion. The simulation package can be used to study different senarios: radiation process, single-defect, defect movement, defect clustering

For each step, the rates of all possible events (F-P genr, vcc jumps, itl jumps) are calculated first, and an event is randomly picked and performed. After that, recombination and annihilation are checked.   

![image](https://github.com/cffrankhuang/W-Re-Kinetic-Monte-Carlo-Simulation/assets/56565245/f43e7fec-bba2-4c90-9067-a3a2ad220849)


### How to use ###

Before use: please install gcc/g++ for the newest version (c++11 or newer)

Set up parameters in "kmc_par.h", and then simply type "make". A executable file named exekmc will be generated. Or compile the source files by "g++ kmc_*.cpp -std=c++11". Run the simulations by using command "./exekmc".

### Physical mechanisms ###
1. Frenkel-pair generation  
    * Two atomic sites are randomly selected, which are then replaced by a vacancy and an interstitial with the type of combination of original atom typess
2. Defect jumping  
    * Vacancy: A vacancy jumps in 1st nearest neghbor. The jump rates are calculated using saddle-point model.
    * Self-interstitial atom(SIA): A SIA is assumed to move in 1-dimension; i.e., it sticks to move in one chosen 1st-nn direction, with certain probability to rotate.
    * mixed interstitial atom: It jumps in 2nd-nn, as an effective jumping mechanism of "bridge-mechanism".


### Program structure ###

* kmc_par.h: store parameters
* kmc_main.cpp: the main file of the program
* kmc_global.h: contains global variables and functions  
    kmc_global.cpp  
* kmc_initial.h: a class used to initialize the simulations  
    kmc_initial.cpp  
* kmc_events.h: a class with all event related variables and functions  
* kmc_events_main.cpp: calculate rates and perform one event  
* kmc_events_genr.cpp: perform F-P generation  
* kmc_events_crcal.cpp: update/calculate vcc creation rate at srf  
* kmc_events_ircal.cpp: calculate itl jumping rates  
* kmc_events_sp.cpp: calculate vcc jumping rates  
* kmc_events_bondecal.cpp: calculate energies using broken bond model  
* kmc_events_recb.cpp: check/perform recombination  

### Developed by ###

Frank Huang. Code adapted from code repository from Marian group at UCLA

