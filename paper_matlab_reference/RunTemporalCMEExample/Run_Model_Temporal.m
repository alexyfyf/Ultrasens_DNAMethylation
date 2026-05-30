function [ModelStruct] = Run_Model_Temporal(Parameters1,Parameters2,DistLength1,DistLength2,tspan)
NCpG=27;
HyperVal=0.8;
HypoVal=0.2;
%the initial distribution is the SS of Parameters1. Then Parameters2
%propagates. Perform for each density value.
ds=[2,3,4,5,6,7,8,9,11,13,17,26,110];
ds=fliplr(ds);

for loopd=1:numel(ds)
    d=ds(loopd)
    CpGPositions=[1:d:NCpG*d];
    Densities=CpGDensities_Function(CpGPositions,50);
    MeanDens(loopd)=mean(Densities);

    %call the function that computes the P(NetMeth) for the coarse-grained
    %approximate model
    [PVecMSM_O,MBins,PVec,Prob_ind,IndCpGp] = CMEModel(NCpG,Parameters1,CpGPositions,DistLength1);
    [PVecMSM_T,MBins] = CMEModel_Temporal(NCpG,Parameters2,CpGPositions,DistLength2,tspan,PVec);

    for loopt=1:numel(tspan)
        PVecMSM=PVecMSM_T(loopt,:);
        %make sure PVec is properly positive and normalized
        PVecMSM(PVecMSM<0)=0;
        PVecMSM=PVecMSM/sum(PVecMSM);
        MethRatio=MBins/NCpG;
        MeanMeth(loopt,loopd)=sum(PVecMSM.*MethRatio);
        HyperInds=find(MethRatio>HyperVal);
        HypoInds=find(MethRatio<HypoVal);
        Hyper(loopt,loopd)=sum(PVecMSM(HyperInds));
        Hypo(loopt,loopd)=sum(PVecMSM(HypoInds));
    end
    %now compute the steady-state for the altered parameters
    loopt=loopt+1;
    [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = CMEModel(NCpG,Parameters2,CpGPositions,DistLength2);
    PVecMSM(PVecMSM<0)=0;
    PVecMSM=PVecMSM/sum(PVecMSM);
    MethRatio=MBins/NCpG;
    MeanMeth(loopt,loopd)=sum(PVecMSM.*MethRatio);
    %HyperInds=find(MethRatio>HyperVal);
    %HypoInds=find(MethRatio<HypoVal);
    Hyper(loopt,loopd)=sum(PVecMSM(HyperInds));
    Hypo(loopt,loopd)=sum(PVecMSM(HypoInds));


end
ModelStruct.Hyper=Hyper;
ModelStruct.Hypo=Hypo;
ModelStruct.Mean=MeanMeth;
ModelStruct.densityvals=MeanDens;
end