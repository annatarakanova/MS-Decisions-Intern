# M&S Decisions summer intern
ABC-SMC Algorithm Settings Testing to Optimize Computations Time for PBPK Model Parameters Optimization

## Introduction 
Physiologically-based pharmacokinetic (PBPK) modeling requires the accurate estimation of unknown drug parameters, a challenge effectively solved by the simulation-based Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm. However, the algorithm's performance critically depends on its initial parameters, which require testing and optimization for the PBPK model of Olaparib. 

## The goal
The main goal was to identify a robust and reproducible set of parameters that balances computational efficiency with the reliable convergence to an accurate posterior parameter distribution.

## The objectives
- To digitize aggregated PK data from the literature and observed dataset preparation.
- To learn basic information about the PBPK model structure and parameters.
- Adaptation of initial R script for particular PBPK model parameter optimization task using ABC-SMC approach. 
- Test the parameters by repeating computations with different ABC-SMC settings and report the results (parameters identifiability, stability of the results, computation time).
**The main question – what are the optimal settings of ABC-SMC approach for a typical PBPK model parameters identification tasks?**

## Process & Methods
**- PBPK modeling**
Mechanistic approach to describe the pharmacokinetics of a substance based on substance-specific properties and mammalian physiology, for which a substantial amount of prior biological information is used for model building. Olaparib is a type of targeted therapy drug called a PARP inhibitor, which pharmacokinetics we were interested in.
**- PK-Sim**
A software tool for whole-body physiologically based pharmacokinetic modeling. This software was used to construct our own project for drug Olaparib and to run the PBPK simulation.
**- ospsuite R package**
A package provides the functionality of loading, manipulating and simulating the simulations created in PK-Sim. It was used to run the simulation by .pkml Olaparib file.
**- Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm**
This iterative algorithm works by proposing parameter sets from a prior distribution, running a PBPK simulation for each set, and accepting only those sets where the simulated output is sufficiently close to the observed clinical data. To serve as gold-standard benchmark, a final reference posterior distribution was generated using a computationally intensive **Rejection from Priors sampling method**.


