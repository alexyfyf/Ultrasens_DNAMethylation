function [p,SSD] = Fitting_MethylationModel(Par0,DL0,DataName)

%Input arguments:
%Par0 - initial parameters for fit
    %Expects 1x14 vector, see CMEMethylationModel for details
%DL0 - initial distance lengthscale parameters for fit
    %Expects 1x12 vector, see CMEMethylationModel for details
%DataName - points to file with data to be fit, expects a structure
    %containing data-derived curve for hypermethylated and hypomethylated
    %CpG fractions as a function of CpG density
    
%Outputs:
%p the vector of fitted parameters
%SSD the SSD between fitted model and data
    
close all
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27; %number of CpGs in the model


DataFile=['Save_' DataName '.mat'];
load(DataFile,'DataStruct')
k1=Par0(1);
k2=Par0(2);
k3=Par0(3);
k13=Par0(13);
k14=Par0(14);

%A subset of the input parameters are selected to be varied in the fitting
p0(1)=Par0(4);
p0(2)=Par0(6);
p0(3)=Par0(8);
p0(4)=Par0(9);
p0(5)=Par0(11);

%add in a little randomness to the initial conditions
p0=p0+1E-3*randn(size(p0));

%initialize the parameters
Parameters=[k1,k2,k3,p0(1),0,p0(2),0,p0(3),p0(4),0,p0(5),0,k13,k14];
DistLength=DL0;
fun=@(p)Fit_Model(p,DataStruct,Parameters)

mns=ones(size(p0))*1E-2;
mxs=[3,25,25,25,25];
options=optimoptions('fmincon','Display','iter');
nonlcon=[];

%run the fit
p = fmincon(fun,p0,[],[],[],[],mns,mxs,nonlcon,options);
Parameters=[k1,k2,k3,p(1),0,p(2),0,p(3),p(4),0,p(5),0,k13,k14]


fn=['FitResult_' DataName];
ffn=[fn '.png']; %figure filename
afn=[fn '.mat']; %array filename

[SSD] = PlotCompare(DataStruct,Parameters,DistLength,ffn);
save(afn,'Parameters','DistLength','SSD')



    function [SSD] = PlotCompare(DataStruct,Parameters,DistLength,fn)
        [ModelStruct] = Run_Model_MorePoints(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);

        dens=ModelStruct.densityvals;
        Hyper=ModelStruct.Hyper;
        Hypo=ModelStruct.Hypo;
        Inter=ones(size(Hyper))-Hyper-Hypo;

        figure(2)
        hyperm=plot(dens,Hyper,'-or')
        hold on
        hypom=plot(dens,Hypo,'-ob')
        plot(dens,Inter,'-ok')

        dens=DataStruct.densityvals;
        Hyper=DataStruct.Hyper;
        Hypo=DataStruct.Hypo;
        Inter=ones(numel(Hyper),1)-Hyper-Hypo;

        hyperd=plot(dens,Hyper,'--r')
        hold on
        hypod=plot(dens,Hypo,'--b')
        plot(dens,Inter,'--k')
        xlabel('Normalized CpG Density')
        ylabel('Fraction of CpGs')
        legend([hyperm,hypom,hyperd,hypod],'Hyper Model','Hypo Model','Hyper Data','Hypo Data')

        text(0.4,0.5,['SSD =' num2str(SSD)])

        print(fn,'-dpng')
    end


    function SSD = Fit_Model(p,DataStruct,Parameters)
        %passes the free fitting parameters into the Parameters vector
        Parameters(4)=p(1);
        Parameters(6)=p(2);
        Parameters(8)=p(3);
        Parameters(9)=p(4);
        Parameters(11)=p(5);
        
        [ModelStruct] = Run_Model(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);
    end


    function SSD = Get_Error(ModelStruct,DataStruct)
        ModelMat=[ModelStruct.Hyper(:),ModelStruct.Hypo(:)];
        %we need to interpolate the experimental data to make it match the
        %model output
        xq=ModelStruct.densityvals;
        Hypervq=interp1(DataStruct.densityvals,DataStruct.Hyper,xq);
        Hypovq=interp1(DataStruct.densityvals,DataStruct.Hypo,xq);
        DataMat=[Hypervq(:),Hypovq(:)];

        Diffs=ModelMat-DataMat;
        SSD=sum(Diffs(:).^2);
    end

    function [ModelStruct] = Run_Model(Parameters,DistLength)
        ds=[2,4,5,6,7,8,9,11,13,17,26,110];
        %ds=[2,4,6,8,11,17,110];
        ds=fliplr(ds);
        
        %loop over different inter-CpG distances to vary the density
        for loopd=1:numel(ds)
            d=ds(loopd);
            CpGPositions=[1:d:NCpG*d];
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec] = CMEMethylationModel(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
    end

    function [ModelStruct] = Run_Model_MorePoints(Parameters,DistLength)
        ds=[2,3,4,5,6,7,8,9,11,13,17,26,110];
        ds=fliplr(ds);

        for loopd=1:numel(ds)
            d=ds(loopd);
            CpGPositions=[1:d:NCpG*d];
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec] = CMEMethylationModel(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
    end

end
