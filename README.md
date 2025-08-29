**run_SALO_par_GMO.jl**: Generates time-series of the phase dynamics for the IEEE 57 test case both with the forcing applied to node ii and without any forcing. The two steps of the algorithm in M. Tyloo, M. Vuffray, A.Y. Lokhov
arXiv preprint arXiv:2310.00458  is applied to the time-series.

**Methods_GMO.jl**: Contains all the functions needed to generate the data and perform step 1 of the algorithm.

**SALO_par_GMO.jl**: Contains step 2 of the algorithm.

**plot_par_IEEE57.jl**: Plots the negative likelihood obtained from the algorithm. 
