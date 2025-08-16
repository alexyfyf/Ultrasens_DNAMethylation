function Plot_ExampleTemporal()

close all

FSS=8; %small font size
FSM=10; %medium font size


NCpG=27; %number of CpGs to simulate
HyperVal=0.8;
HypoVal=0.2;


%load the baseline parameters (WT fit) and plot 
Name='Fit_HUES8WT_CpGsOnly_Chr1_58224'; 
load(Name,'Parameters','DistLength')
[ModelStruct] = Run_Model(Parameters,DistLength); %run the steady-state CME model according to loaded parameters
[ED50,nslope,KFD,nFD,KFM,nFM] = FitSwitch_MultMethods(ModelStruct); %fit the methylation switch to the model output
BaseVals=[ED50,nslope,KFD,nFD,KFM,nFM];

%call the function that runs the temporal model and plots it
PlotExample(Parameters,DistLength);

figure(1)
 set(gcf,'Units','inches');
        outpos=get(gcf,'Position');
        outpos(3)=2.75;%7.75%7.5;
        outpos(4)=2;
        set(gcf,'Position',outpos)
        set(gca,'FontSize',FSS)
        %axis square
        xlabel('Local CpG Density')
        ylabel('Mean CpG Methylation')
        set(gca,'FontSize',FSM)
        xlim([0 0.5])

fontname("Times")

print -dpng ParameterSweepsExample
print -dtiff ParameterSweepExample

figure(2)
fontname("Times")
print -dpng Inset
print -dtiff Inset

    function PlotExample(Parameters,DistLength)
 
        %change the parameters somehow
        NewPar=Parameters;
        NewDist=DistLength;
        NewPar(3)=Parameters(3)*0.88;
        %set the time-points for the temporal model
        tspan=[[0:20:60],100,200];

        %call the function that runs the temporal CME model. It is
        %initialized with the steady-state distribution corresponding to
        %[Parameters,DistLength]. Then it evolves the CME model in time
        %according to a modified set of parameters in [NewPar,NewDist]
        [ModelStruct] = Run_Model_Temporal(Parameters,NewPar,DistLength,NewDist,tspan);
        ModelStructTemp=ModelStruct;
  
        ModelStruct=ModelStructTemp;
        HyperA=ModelStruct.Hyper;
        numt=size(HyperA,1);
       
        figure(1)
   
        hold on
        mygry=[0.6 0.6 0.6];
        %ploy the temporal model output
        for tind=1:numt-1%numel(keepinds)
            %tind=keepinds(loopt);
            dens=ModelStruct.densityvals;
            Mean1=ModelStruct.Mean(tind,:);
            Hyper1=ModelStruct.Hyper(tind,:);
            Hypo1=ModelStruct.Hypo(tind,:);
            %Clr=ColorsS(cinds(loopt),:);
            yy=smooth(Mean1);
            if tind==1
                plot(dens,yy,'-k','LineWidth',2)
            else
                plot(dens,yy,'-','LineWidth',2,'Color',mygry)
            end
        end
        %fit the methylation switch at the various timepoints
        [AllVals,numt] = FitSwitch_Temporal(ModelStruct);
        KMList=AllVals(:,5)-BaseVals(5);%ED50
        nMList=AllVals(:,6)-BaseVals(6);%nslope

        figure(2)
        bubblechart(KMList(1:end-1),nMList(1:end-1),1,'MarkerEdgeColor',mygry,...
            'MarkerFaceColor',mygry,'MarkerFaceAlpha',0.2, ...
            'MarkerEdgeAlpha',0.6)

        hold on
        bubblechart(KMList(1),nMList(1),1,'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.2, ...
            'MarkerEdgeAlpha',0.6)
        %plot(KMList,nMList,'-k','LineWidth',2)
        figure(2)
        bubblesize(gca,[5 5.1])


        figure(2)
        axis square


        grid on
        xlim([-0.1 0.1])
        ylim([-2 2])
        set(gcf,'Units','inches');
        outpos=get(gcf,'Position');
        outpos(3)=1.1;%7.75%7.5;
        outpos(4)=1.1;
        set(gcf,'Position',outpos)
        set(gca,'FontSize',FSS)
        xlabel('K*','FontSize',FSM)
        ylabel('n*','FontSize',FSM)

    end


    function [ED50,nslope,KFD,nFD,KFM,nFM] = FitSwitch_MultMethods(ModelStruct)
        Hypo=ModelStruct.Hypo;
        Hyper=ModelStruct.Hyper;
        Mean=ModelStruct.Mean;
        xax=ModelStruct.densityvals;
        Diff=Hyper-Hypo;
        CDiff=Diff+1;
        Dirct=0;

        [ED50,nslope,flag] = GetED50(Hypo,Hyper,xax);
        [nFD,KFD,HillModel,Err] = FitHill2Param(xax,CDiff,Dirct);
        [nFM,KFM,HillModel,Err] = FitHill2Param(xax,Mean,Dirct);

    end

    function [AllVals,numt] = FitSwitch_Temporal(ModelStruct)
        HyperA=ModelStruct.Hyper;
        numt=size(HyperA,1);
        AllVals=zeros(numt,6);
        for tind=1:numt
            dens=ModelStruct.densityvals;
            Mean1=ModelStruct.Mean(tind,:);
            Hyper1=ModelStruct.Hyper(tind,:);
            Hypo1=ModelStruct.Hypo(tind,:);
            Diff=Hyper1-Hypo1;
            CDiff=Diff+1;
            Dirct=0;

            [ED50,nslope,flag] = GetED50(Hypo1,Hyper1,dens);
            [nFD,KFD,HillModel,Err] = FitHill2Param(dens,CDiff,Dirct);
            [nFM,KFM,HillModel,Err] = FitHill2Param(dens,Mean1,Dirct);
            AllVals(tind,:)=[ED50,nslope,KFD,nFD,KFM,nFM];
        end

    end


    function [ED50,slope,flag] = GetED50(Hypo,Hyper,xax)

        % % %no interpolation method:
        % AAA=Hypo-Hyper;
        % getinds=find(AAA>-0.3 & AAA<0.3);
        % ED50=interp1(AAA(getinds),xax(getinds),0);
        % flag=isnan(ED50);
        % slopesfine=gradient(AAA)./gradient(xax);
        % slope=interp1(xax,slopesfine,ED50);
        %
        % %another fit with a trendline (and interpolation is necessary)
        % %Seems to be almost identical to point-estimate of slope, above
        % finex=linspace(min(xax),max(xax),300);
        % Hypofine=interp1(xax,Hypo,finex);
        % Hyperfine=interp1(xax,Hyper,finex);
        % AAA = (Hypofine-Hyperfine);
        % getinds=find(AAA>-0.3 & AAA<0.3);
        % p=polyfit(finex(getinds),AAA(getinds),1);
        % slope=p(1);


        %%%%%with interpolation method:
        finex=linspace(min(xax),max(xax),300);
        Hypofine=interp1(xax,Hypo,finex);
        Hyperfine=interp1(xax,Hyper,finex);
        AAA = (Hypofine-Hyperfine);

        %first method, using interp1
        getinds=find(AAA>-0.2 & AAA<0.2);
        if numel(getinds)>3
            ED50=interp1(AAA(getinds),finex(getinds),0);
            p=polyfit(finex(getinds),AAA(getinds),1);
            slope=p(1);
            flag=0;
        else
            ED50=0;
            slope=0;
            flag=1;
        end
        %
        % %another method to calculate ED50 (doesn't seem to matter which)
        % %[v,inn]=min(abs(AAA));
        % %ED50=finex(inn);
        % %if v>0.05
        % %    flag=1; %flag if it doesn't cross
        % %else
        % %    flag=0;
        % %end
        %
        % %another method to calculate slope (after interpolation)
        % slopesfine=gradient(AAA)./gradient(finex);
        % slope=interp1(finex,slopesfine,ED50);

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
        ModelStruct.Mean=MeanM;
    end
end
