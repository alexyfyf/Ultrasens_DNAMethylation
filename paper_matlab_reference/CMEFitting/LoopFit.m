clear
close all
NLoops=3; %number of times to run the fit
Nt=24; %cell cycle time in hours
k13=1/Nt;
k14=1/(2*Nt);
for ii=1:NLoops
    ii

%initialize the parameters (see Fit_CME_Methylation for parameter
%definitions)
Parameters=[0.0999, 0.1000, 1, 0.1431, 0, 3.5747, 0, 0.7281, ...
     0.0639, 0, 24.7399, 0, k13, k14];
DistLength= [0, 0, 0, 0, 33.6308, 33.6308, 33.6308, 33.6308, ...
    43.3724, 43.3724, 43.3724, 43.3724];

DataName='HUES8WT_CpGsOnly_Chr1';

     [p,SSD] = Fit_CME_Methylation_PS(Parameters,DistLength,DataName)

end