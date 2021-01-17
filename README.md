# power-aware-task-scheduling

This is part of a project under the supervision of [Prof. Dr. Abhishek Mishra](https://www.bits-pilani.ac.in/pilani/abhishekmisra/profile). The objective is to analyse and come up with algorithms that schedule tasks in certain parallel computing environments given an Energy constraint.

The project is divided into two parts.  
1. The first is an implementation of the Energy contraint task scheduling algorithms that minimize Schedule length and maximize Reliability, the algorithms of which have been discussed in the [following paper](https://ieeexplore.ieee.org/document/8936469). The results produced by the implementations are exactly the same as mentioned  in the paper. **Look at ISAECC, ISAECCStar and MRDECC for the implementations.**  

2. The second part is a proposed Monte Carlo algorithm that schedules tasks given an energy constraint in an environment similar to the one described in the [following paper](https://www.sciencedirect.com/science/article/pii/S0307904X13006355). Along with the implementation of this Monte Carlo Algorithm, an optimum(exponential) algorithm with certain optimizations for faster running time has also been implemented for comparison with the first algorithm. **Look at MonteCarloAlgorithm and OptimumAlgorithm for the implementations.**  
