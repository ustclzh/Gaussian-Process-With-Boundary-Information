These codes provide MATLAB implementation of examples in Section 4 of the paper "Improving Gaussian Process Emulators with Boundary Information ". 

In particular,
Implementation of example in Section 4.1: Folder named 'example4.1';
Implementation of example in Section 4.2: Folder named 'example4.2'.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To obtain results in Section 4.1, open the folder titled "example4.1" and then run the .m file 'main.m' ; 
To obtain results in Section 4.2, open the folder titled "example4.2" and then run the .m files：
       run 'one_dim_boundary_inf.m', you will obtain figure 4 in section 4.
       run 'one_dim_boundary_r.m', you will obtain figure 5 in section 4.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
After MATLAB has finished with the computations, the results are stoared in the corresponding root folder:
For example 4.1: after runing 'main.m', the result is stored by .mat file 'RESULT.mat' in folder named 'example4.1': 'result.mat' includes the results shown in table 1. In particular, 
res.summary: Sample mean (in parentheses) of the MAE, RMSE, ALPI, and Coverage given by the standard GP, FBGP, BdryGP, BMGP, and GGP emulators for 60 pairs of design and test sets
res.se: Standard error (in parentheses) of the MAE, RMSE, ALPI, and Coverage given by the standard GP, FBGP, BdryGP, BMGP, and GGP emulators for 60 pairs of design and test sets

For example 4.2: after runing 'one_dim_boundary_inf.m', the results are stored by .tif file 'result_boundary_inf.tif' in folder named 'example4.2': 'result_boundary_inf.tif' is the figure 4;
	          after runing 'one_dim_boundary_r.m', the results are stored by .tif file'result_boundary_r.tif' in folder named 'example4.2': 'result_boundary_r.tif' is the figure 5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



