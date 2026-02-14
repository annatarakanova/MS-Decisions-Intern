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
- PBPK modeling

Mechanistic approach to describe the pharmacokinetics of a substance based on substance-specific properties and mammalian physiology, for which a substantial amount of prior biological information is used for model building. Olaparib is a type of targeted therapy drug called a PARP inhibitor, which pharmacokinetics we were interested in.

- PK-Sim

A software tool for whole-body physiologically based pharmacokinetic modeling. This software was used to construct our own project for drug Olaparib and to run the PBPK simulation.

- ospsuite R package

A package provides the functionality of loading, manipulating and simulating the simulations created in PK-Sim. It was used to run the simulation by .pkml Olaparib file.

- Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) algorithm

This iterative algorithm works by proposing parameter sets from a prior distribution, running a PBPK simulation for each set, and accepting only those sets where the simulated output is sufficiently close to the observed clinical data. To serve as gold-standard benchmark, a final reference posterior distribution was generated using a computationally intensive **Rejection from Priors sampling method**.

## Results
- The aggregated clinical PK data were digitized from published articles and observed datasets were prepared – Olaparib 300 mg SD (single dose) from Plummer et al. (2015) and Olaparib 400 mg BID (twice a day) from Gao et al. (2023). 
- The similarity metric was created to show the similarity/difference between posterior parameter distributions that were obtained by Rejection from Priors sampling method and ABC-SMC approach. 
- 10 parameters of ABC-SMC algorithm were tested: number of iterations (Niterations) = 4, 6, 9; start epsilon (start_epsilon) = 0.3, 0.5, 1.0; number of solutions after the 1st iteration (Nposterior) = 20, 50, 100; step size while finding new solutions (ScCov) = 0.25, 0.5, 1.0; another observed dataset (400_Gao_2023) with base metrics/start_epsilon = 0.3.
- Graphs of posterior parameter distribution comparison were obtained for each of the computation with different sets of parameters. 
- Reproducibility of results based on base metrics (id = 0) and start_epsilon = 0.3 (id = 2a) were tested.

## Conclusion
In this work, we systematically optimized ABC-SMC algorithm parameters to best match a reference posterior distribution generated via Rejection Sampling. An optimal configuration was identified that provides the best balance of accuracy, speed and reproducibility: a strict starting epsilon (~0.3-0.5), number of solutions after the 1st iteration  (Nposterior) of about 50, a small covariance coefficient (ScCov) of 0.25 and 4-6 iterations. The study also revealed that model parameters were not fully transferable between different clinical datasets, indicating potential model or data-specific limitations. Moreover, the results of reproducibility test indicate that the further research is needed. 

Future research should analyze the initial parameter sets from the first iteration to identify features that correlate with the final stability of the algorithm. Another direction is to explore the impact of sampling parameters in log-space on the overall efficiency and convergence of the optimization process. Finally, a comprehensive statistical analysis of the reproducibility for the most promising scenarios should be conducted to rigorously validate that their high performance is consistent and not due to chance.




