%% CGI-Level Analysis, Human Cell Types
clear
clc
close all
% Load data for Island WGBS Annotations. Each CGI is treated as 1 data point.
DataNa = [{['IslandLvl_agg']},{'IslandLvl_agg_HUES8'},{'IslandLvl_agg_IMR90'}]; % HUES64, HUES8, IMR90 Data loading for loop
TitleNA = [{'HUES64'},{'HUES8'},{'IMR90'}];

% Analysis Parameters
Hypocutoff=0.2; %cutoff below which island is "hypo"
Hypercutoff=1-Hypocutoff;
%Hypercutoff=0.8; %cutoff above which island is "hyper"
fonttt = 24;
lenLim = [205 4500];
linsz=2;

Bivarxticks=[205 2000 4000];
CGIlims=[205 750 1500];

% Initialize Subplot axes
figure(1000)
for k= 1:(2*numel(DataNa))
	ax1(k)=subplot(2,numel(DataNa),k);
end

%% Bivariate Histogram Alg.
% As a first pass, we will be comparing/looking at the dependence of methyl
% state on CGI Length

for j= 1:(numel(DataNa))
	M4 = csvread(DataNa{j},1,1); % Load Data for the loop
	
	x = M4(:,4); % CGI Length in bp
	x1 = M4(:,3); % Number of CpGs in associated CGI
	y = M4(:,2); % Mean WGBS Value for each CGI
	N = numel(x); % This should give us the total number of CGIs/Datapoints
	
	xgridS=linspace(0,max(x),500);
	xgrid1=linspace(0,max(x1),500);
	ygrid=linspace(0,max(y),100);
	
	NumBinsX=numel(xgridS);
	NumBinsX1=numel(xgrid1);
	NumBinsY=numel(ygrid);
	
	%initialize the Bivariate array
	Bivariate=zeros(NumBinsX,NumBinsY);
	Bivariate_c=zeros(NumBinsX1,NumBinsY);
	
	%loop through and count members in each gridpoint
	for ii=1:N
    	indsx=find(x(ii)<=xgridS);
    	indsy=find(y(ii)<=ygrid);
    	Bivariate(indsx(1),indsy(1))=Bivariate(indsx(1),indsy(1))+1; %increment the array in this location
	end
	indsx=[];
	indsy=[];
	for k=1:N
    	indsx=find(x1(k)<=xgrid1);
    	indsy=find(y(k)<=ygrid);
    	Bivariate_c(indsx(1),indsy(1))=Bivariate_c(indsx(1),indsy(1))+1; %increment the array in this location
	end
	
	% Classification Scheme for Hyper-/Hypomethylated Islands:
	
	%To bin the lengths, use bins based on percentiles. This will give approx
	%same number of islands per bin. (I'm not quite sure why it's not EXACTLY
	%the same number of islands per bin.
	ygrid3=0:0.05:1;
	Nbins=50; %Number of bins in length
	Percs=linspace(0,100,Nbins); %these will be the requested percentiles between 0 and 100
	xgrid=prctile(x,Percs); %these are the Bin edges for CGIlength
	BivariateDist=zeros(numel(xgrid),numel(ygrid3)); %initialize the array
	
	%build the bivariate heat map w/Percentile Bin Edges
	indsx=[];
	indsy=[];
	for ii=1:N
    	indsx=find(x(ii)<xgrid);
    	indsy=find(y(ii)<ygrid3);
    	if indsx
    	else
        	indsx=size(BivariateDist,1);
    	end
    	if indsy
    	else
        	indsy=size(BivariateDist,2);
    	end
    	BivariateDist(indsx(1),indsy(1))=BivariateDist(indsx(1),indsy(1))+1;
	end
	
	% Using a Cutoff/Threshold value, assign islands to be in the hyper-,
	% hypo-, or intermediate methylation state.
	
	Hypobins=find(ygrid3<Hypocutoff);
	Hyperbins=find(ygrid3>Hypercutoff);
	Interbins=setdiff(1:numel(ygrid3),[Hyperbins,Hypobins]);
	
	%lump the WGBS bins together into the three classes
	%(cut off the last CGIlegnth bin for graphing purposes)
	Hypoarray=sum(BivariateDist(2:end-1,Hypobins),2);
	Interarray=sum(BivariateDist(2:end-1,Interbins),2);
	Hyperarray=sum(BivariateDist(2:end-1,Hyperbins),2);
	xax=xgrid(2:end-1);
	%% Visualization
	Bivariate2=Bivariate/numel(x); % Normalize CGI Length bivariate data to number of CGIs
	Bivariate2_c=Bivariate_c/numel(x); % Normalize data CpG Number bivariate data to number of CGIs
	
	subplot(ax1(j))
	s1=pcolor(xgridS,ygrid,[-log(Bivariate2)]');
	colormap(flipud(parula))
	yline(Hypocutoff,'--g','LineWidth',linsz)
	yline(Hypercutoff,'--r','LineWidth',linsz)
	xlabel('CGI Length (bp)')
	ylabel('Methyl. Frac.')
	xlim([lenLim])
	xticks(Bivarxticks)
	title(TitleNA{j})
	ax=gca;
	ax.FontSize=fonttt-4;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];
	
	subplot(ax1(j+3))
	plot(xax,Hypoarray,'-o','LineWidth',linsz,'Color','b')
	hold on
	plot(xax,Hyperarray,'-o','LineWidth',linsz,'Color','r')
	plot(xax,Interarray,'-o','LineWidth',linsz,'Color','k')
	xlabel('CGI Length')
	ylabel('No. of Islands')
	%xlim([min(xax) max(xax)])
	xlim([205 1500])
	ylim([0 100])
	xticks(CGIlims)
	yticks([0 250 500])
	ax=gca;
	ax.FontSize=fonttt-4;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];

	%Calculating Hypo-Hyper Cross over point
	AAA = (Hypoarray-Hyperarray);
	[v,inn]=min(abs(AAA));
	title(['CGI-Cross = ' num2str(xax(inn)) 'bp'],'FontSize',(fonttt-4))

end
subplot(ax1(numel(DataNa)))
h=colorbar;
ylabel(h,'-log(Probability)')
yline(Hypocutoff,'--g','Hypomethyl.','LineWidth',linsz,'FontSize',(fonttt-12))
yline(Hypercutoff,'--r','Hypermethyl.','LineWidth',linsz,'FontSize',(fonttt-12))

subplot(ax1(2*numel(DataNa)))
legend(['Hypo < ' num2str(Hypocutoff)],['Hyper > ' num2str(Hypercutoff)],'Inter','Location','best','FontSize',(fonttt-4))