clear
close all
%this file loads two different data files and plots data about Delta
%Methylation for individual CGIs between the two cell types
CompFiles = {'IslandLvl_agg_HUES64WT','IslandLvl_agg_HUES64_DKO_Early','IslandLvl_agg_HUES64_DKO_Late'};
CompNAs={'WT','DKO (Early)','DKO (Late)'};

%% Explanation
% This file compares the HUES64WT to its DKO at early/timepoints in a
% hardcoded way.
indx=[1,2; 1,3; 2,3];
subLoc=[1,2,3.5];
for zz=1:3

%% File Loading
	%set the filenames and Cell Names (make sure they're consistent!)
	filename=CompFiles{indx(zz,1)};
	CellName1=CompNAs{indx(zz,1)};
	Data1=readtable(filename);
	
	filename=CompFiles{indx(zz,2)};
	CellName2=CompNAs{indx(zz,2)};
	Data2=readtable(filename);
	
	numisles=numel(Data1.CGIno);
	
	% Begin Analysis
	
	ct=0;
	%need to find all the islands that are common between the two datasets
	%careful! It does this by checking length and CGIno, but it would be better
	%if the location in the genome were used as a unique identifier for CGIs
	%across cell types (but this way probably works, though it's possible some
	%were misidentified across the files because two different CGIs could have
	%same length and number). 
	
	
	for ii=1:numisles
    	%ii=1; %this is the index in the first dataset
    	CpGNum1=Data1.CpGNum(ii);
    	CGIlen1=Data1.CGIlen(ii);
	
    	Numinds=find(Data2.CpGNum==CpGNum1);
    	Leninds=find(Data2.CGIlen(Numinds)==CGIlen1);
    	Inds=Numinds(Leninds);
	
    	%if there are more than one island in the other data set with same
    	%#/length, choose the one thats closer in the file (closer in terms of
    	%the CGIno)
	
    	diffinds=Data1.CGIno(ii)-Data2.CGIno(Inds);
    	[val,ind]=min(abs(diffinds));
	
    	if ind
        	ct=ct+1;
	
        	jj=Inds(ind); %this is the index in the second dataset
        	DeltaMeth(ct)=Data1.WGBS(ii)-Data2.WGBS(jj);
        	IndsList(ct,:)=[ii,jj];
    	end
	end
	CommonInds=IndsList(:,1);
	
	numisles=ct; %revise the number CGIs to just the number in common between datasets
	Lengths=Data1.CGIlen(CommonInds); %keep the common lengths
	%bin by percentiles
	Percentiles=[0:2:100];
	for pr=1:numel(Percentiles)
    	LenBin(pr)=prctile(Lengths,Percentiles(pr));
	end
	
	[N,LenBin,bin]=histcounts(Lengths,LenBin);
	minb=min(bin);
	maxb=max(bin);
	MethDiffBins=[-1,-0.6,-0.2,0.2,0.6,1];
	for bb=minb:maxb
    	getinds=find(bin==bb);
    	Meanbin(bb)=mean(Lengths(getinds)); %the mean length in this bin
    	MeanDiff(bb)=mean(DeltaMeth(getinds)); %the mean change of CGI meth in this bin
    	MedDiff(bb)=median(DeltaMeth(getinds)); %the mean change of CGI meth in this bin
    	Prc25(bb)=prctile(DeltaMeth(getinds),25);
    	Prc75(bb)=prctile(DeltaMeth(getinds),75);
    	%Prc95(bb)=prctile(DeltaMeth(getinds),95);
    	NumHere=numel(getinds);
    	
	%bin the CGIs into difference categories
	[N,MethDiffBins,delbin]=histcounts(DeltaMeth(getinds),MethDiffBins);
	Classes(bb,:)=N/NumHere*100;
	
	end
	%% Kojo Edit Jan 30: Fig2, Except we ignore the ones with minimal change.
	figure(3000)
	%subplot(2,2,subLoc(zz))
	subplot(1,3,zz)
	[at,bt]= size(Classes);
	Classes2=zeros(at,(bt-1));
	Classes2(:,1:2)=Classes(:,1:2);
	Classes2(:,3:4)=Classes(:,4:5);
	semilogx(Meanbin,Classes2,'-','LineWidth',3)
	
	newcolors=[1 0.3 0.3;
    	0.6 0 0;
    	0 0 0.6;
    	0.3 0.3 1];
	colororder(newcolors)
	
	xlim([200 max(Meanbin)])
	axis square
	xlabel('CGI Length (bp)')
	ylabel('% of CGIs')
	title([CellName1 ' vs. ' CellName2])
	set(gca,'FontSize',14)
	grid on
	xticks([200 500 1000 2000])
	yticks([0 10 20])
end

% NOTE: To create the figure, the legend must be manually moved in the GUI.