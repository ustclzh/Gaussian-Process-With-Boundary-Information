# Gaussian-Process-With-Boundary-Information


These codes provide MATLAB implementation of examples in Section 4 of the paper "Improving Gaussian Process Emulators with Boundary Information ". 

In particular,
Implementation of example in Section 4.1: Folder named 'example4.1';
Implementation of example in Section 4.2: Folder named 'example4.2'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

To obtain results in Section 4.1, open the folder titled "example4.1" and then run the .m file 'main.m' ; 
To obtain results in Section 4.2, open the folder titled "example4.2" and then run the .m files 'one_dim_boundary_inf.m' and 'one_dim_boundary_r.m' respectively. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

After MATLAB has finished with the computations, the results are stoared in the corresponding root folder:
For example 4.1: after runing 'main.m', the result is stored by .mat file 'RESULT.mat' in folder named 'example4.1': 'RESULT.mat' includes the results shown in table 1. 
For example 4.2: after runing 'one_dim_boundary_inf.m', the results are stored by .tif file 'result_boundary_inf.tif' in folder named 'example4.2': 'result_boundary_inf.tif' is the figure 4;
	          after runing 'one_dim_boundary_r.m', the results are stored by .tif file'result_boundary_r.tif' in folder named 'example4.2': 'result_boundary_r.tif' is the figure 5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

When MATLAB finishes, the variables below are important. 
For example 4.1: the variable 'res' is a summary of result same as in Table 1 of paper. 
