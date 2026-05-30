function [p,SSD] = Fit_CME_Methylation_PS(Par0,DL0,DataName)

%Par0 - initial parameters for fit
%[1x14] vector. Indexing of reactions:
%u+,h-,h+,m-,h+h,h+m,u+h,u+m,h-u,h-h,m-u,m-h
%13=m-(passive, replication), 14=h- (passive, replication)
%DL0 - initial distance lengthscale parameters for fit
%[1x12] vector. Indexing: lengthscales of reactions 1-12. (Elements 1-4 are
%not used, since these are the standard reactions)
%DataName - points to right dataset to be fit
close all
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27; %number of CpGs in the model


DataFile=['Save_' DataName '.mat'];
load(DataFile,'DataStruct')

%initializes the parameters to be fit
p0(1)=Par0(1);
p0(2)=Par0(2);
p0(3)=Par0(4);
p0(4)=Par0(6);
p0(5)=Par0(8);
p0(6)=Par0(9);
p0(7)=Par0(11);
p0(8)=DL0(5);
p0(9)=DL0(9);

%initialize the parameters
Parameters=Par0;
DistLength=DL0;
fun=@(p)Fit_Model(p,DataStruct,Parameters,DistLength)%,CpGPosArray)

%set the minimum and maximum range for fitted parameters
mns=ones(size(p0))*1E-5;
mxs=ones(size(p0))*25;
mxs(1:3)=1;
mxs(8:9)=100;
nvars=numel(p0);

if exist('savex.mat')
    load('savex.mat','savex')
else
    savex=[];
end
%options=optimoptions('fmincon','Display','iter');
options = optimoptions('particleswarm','Display','iter','HybridFcn',@fmincon, ...
    'InitialSwarm',savex,'OutputFcn',@outfun);
%nonlcon=[];
%run the fit
%p = fmincon(fun,p0,[],[],[],[],mns,mxs,nonlcon,options);
p = particleswarm(fun,nvars,mns,mxs,options);
p
[Parameters,DistLength]=UpdateParameters(p,Parameters,DistLength)

tag=randi(100000); %generate a random integer tag. This will be a unique identifier for the output
%from the fitting
fn=['Fit_' DataName '_' num2str(tag)];
ffn=[fn '.png']; %figure filename
afn=[fn '.mat']; %array filename

[SSD] = PlotCompare(DataStruct,Parameters,DistLength,ffn);%,CpGPosArray);
save(afn,'Parameters','DistLength','SSD')

    function [SSD] = PlotCompare(DataStruct,Parameters,DistLength,fn)%,CpGPosArray)
        %[ModelStruct] = Run_Model_CallPoints(Parameters,DistLength,CpGPosArray);
        [ModelStruct] = Run_Model_MorePoints(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);

        dens=ModelStruct.densityvals;
        Hyper=ModelStruct.Hyper;
        Hypo=ModelStruct.Hypo;
        Inter=ones(size(Hyper))-Hyper-Hypo;

        figure(2)
        plot(dens,Hyper,'-or')
        hold on
        plot(dens,Hypo,'-ob')
        plot(dens,Inter,'-ok')

        dens=DataStruct.densityvals;
        Hyper=DataStruct.Hyper;
        Hypo=DataStruct.Hypo;
        Inter=ones(numel(Hyper),1)-Hyper-Hypo;

        plot(dens,Hyper,'--r')
        hold on
        plot(dens,Hypo,'--b')
        plot(dens,Inter,'--k')

        text(0.4,0.5,['SSD =' num2str(SSD)])

        print(fn,'-dpng')
    end


    function SSD = Fit_Model(p,DataStruct,Parameters,DistLength)%,CpGPosArray)
        [Parameters,DistLength]=UpdateParameters(p,Parameters,DistLength);
        %[ModelStruct] = Run_Model_CallPoints(Parameters,DistLength,CpGPosArray);
        [ModelStruct] = Run_Model(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);
    end


    function SSD = Get_Error(ModelStruct,DataStruct)
        ModelMat=[ModelStruct.Hyper(:),ModelStruct.Hypo(:),ModelStruct.MeanMeth(:)];
        %we need to interpolate the experimental data to make it match the
        %model output
        xq=ModelStruct.densityvals;
        Hypervq=interp1(DataStruct.densityvals,DataStruct.Hyper,xq);
        Hypovq=interp1(DataStruct.densityvals,DataStruct.Hypo,xq);
        Meanvq=interp1(DataStruct.densityvals,DataStruct.MeanMeth,xq);
        DataMat=[Hypervq(:),Hypovq(:),Meanvq(:)];

        Diffs=ModelMat-DataMat;
        SSD=sum(Diffs(:).^2);
    end

   function [ModelStruct] = Run_Model_CallPoints(Parameters,DistLength,CpGPosArray)
        %load('CpGPosArray.mat','CpGPosArray');
        numpoints=size(CpGPosArray,1);
        NCpG=size(CpGPosArray,2);

        for loopd=1:numpoints

            CpGPositions=CpGPosArray(loopd,:);
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = CMEModel(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
            MeanM(loopd)=sum(PVecMSM.*MethRatio);
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
        ModelStruct.MeanMeth=MeanM;
    end

    function [ModelStruct] = Run_Model(Parameters,DistLength)
        ds=[2,4,5,6,7,8,9,11,13,17,26,110];
        %ds=[2,4,6,8,11,17,110];
        ds=fliplr(ds);

        for loopd=1:numel(ds)
            d=ds(loopd);
            CpGPositions=[1:d:NCpG*d];
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = CMEModel(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
             MeanM(loopd)=sum(PVecMSM.*MethRatio);
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
         ModelStruct.MeanMeth=MeanM;
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
            [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = CMEModel(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
             MeanM(loopd)=sum(PVecMSM.*MethRatio);
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
         ModelStruct.MeanMeth=MeanM;
    end

    
    function [Parameters,DistLength] = UpdateParameters(p,Parameters,DistLength)
        Parameters(1)=p(1);
        Parameters(2)=p(2);
        Parameters(4)=p(3);
        Parameters(6)=p(4);
        Parameters(8)=p(5);
        Parameters(9)=p(6);
        Parameters(11)=p(7);
        DistLength(5:8)=p(8);
        DistLength(9:12)=p(9);
    end

    function stop = outfun(optimValues,state)
        %this function stops the optimization after T seconds
        % T=400;
        % stop = toc>T;
        % if stop
        %     toc
        % end
        stop=0;
        switch state
            case 'iter'
                itr=optimValues.iteration;
                if mod(itr,10)==0
                    savex=optimValues.bestx;
                    save savex savex
                end

            case 'done'
        end
    end

end
