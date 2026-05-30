%points to data file
DataName='HUES8WT_CpGsOnly_Chr1';

%initialize some parameters
Nt=24; %cell cycle length in hours
k13=1/Nt;
k14=1/(2*Nt);
Par0=[0.1,0.1,4,0.14,0,3.5,0,0.73,0.06,0,25,0,k13,k14];
DL0=[0,0,0,0,33.63,33.63,33.63,33.63,43.38,43.38,43.38,43.38];

%Call the fitting code
[p,SSD] = Fitting_MethylationModel(Par0,DL0,DataName)
