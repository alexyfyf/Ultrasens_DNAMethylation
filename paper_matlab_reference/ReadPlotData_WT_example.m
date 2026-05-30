clear all
close all

%% Data Loading
% KB Edit: Loading Data Files in an Automated way
%Get directory folder info
% Other folders just add "_Wind##"
FoldExt='./Proc_CSVs/'; % Name of folder with files from Python/Pandas; this is initially called "MATLAB_Proc_CSVs.." then the filed are moved over after converting them from .bed to .csv
F = dir('Proc_CSVs');

%Checking Out the Window Sizes
% FoldExt='./Proc_CSVs_Wind100/'; % Name of folder with files from Python/Pandas; this is initially called "MATLAB_Proc_CSVs.." then the filed are moved over after converting them from .bed to .csv
% F = dir('Proc_CSVs_Wind100');

% FoldExt='./Proc_CSVs_Wind200/'; % Name of folder with files from Python/Pandas; this is initially called "MATLAB_Proc_CSVs.." then the filed are moved over after converting them from .bed to .csv
% F = dir('Proc_CSVs_Wind200');

% FoldExt='./Proc_CSVs_Wind1000/'; % Name of folder with files from Python/Pandas; this is initially called "MATLAB_Proc_CSVs.." then the filed are moved over after converting them from .bed to .csv
% F = dir('Proc_CSVs_Wind1000');
%names={F(:).name};

names{1}='HUES64WT_CpGsOnly_Chr1';
names{2}='HUES8WT_CpGsOnly_Chr1';
names{3}='IMR90WT_CpGsOnly_Chr1';

%Suffix to cut off
cut='_CpGsOnly_Chr1';

nameN=cell(numel(names),1);
for i=1:numel(names)
	a=names{i};
	nameN{i}=erase(a,cut);
end
currdir=pwd;

% Analysis/Plotting Parameters
Hypocutoff=0.2; %cutoff below which CpG is "hypo"
Hypercutoff=1-Hypocutoff;
fonttt = 24;
linsz=2;

%Data=csvread(DataFile,1,1);

DataNa=names;
TitleNA=nameN;
QuanCr=zeros(numel(names),1);
numNN=length(names);
fs=20;
CGIlentick=[0 0.2 0.4];
%% Loop Initialization
for j= 1:(numel(DataNa))
	Data = csvread([FoldExt DataNa{j}],1,1); % Load Data for the loop
	%% Begin Analysis
	x = Data(:,4); % Calculated CpG density
	y = Data(:,3); % WGBS Value for each CpG
	N = numel(x); % This should give us the total number of CpGs/Datapoints
	
	%just in case the methylation values are NOT normalized to 1, for some
	%reason:
	y=y/max(y);
	densityvals=unique(x);
	
	%set up the bins in x and y
	ygrid=linspace(0,max(y),50); %methylation bins
	xgridS=densityvals; %we just use the raw density values--they're already discretized by nature
	
	NumBinsX=numel(xgridS);
	NumBinsY=numel(ygrid);
	
	%initialize the Bivariate array
	Bivariate=zeros(NumBinsX,NumBinsY);
	
	%loop through and count members in each gridpoint to construct
	%bivariate heatmap
	for ii=1:N
    	indsx=find(x(ii)<=xgridS);
    	indsy=find(y(ii)<=ygrid);
    	Bivariate(indsx(1),indsy(1))=Bivariate(indsx(1),indsy(1))+1; %increment the array in this location
	end
	
	%% Visualization
	Bivariate2=Bivariate/numel(x); % Normalize bivariate data to plot density/probability, instead of counts
	
	%plot the bivariate distribution
	figure(1000)
	%subplot(ax1(j))
	subplot(2,numNN,j)
	s1=pcolor(xgridS,ygrid,[-log(Bivariate2)]');
	colormap(flipud(parula))
	yline(Hypocutoff,'--g','LineWidth',linsz)
	yline(Hypercutoff,'--r','LineWidth',linsz)
	%xlabel('CGI Length (bp)')
	ylabel('Methyl. Frac.')
	ax=gca;
	ax.FontSize=fs;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];
	yticks([0 0.5 1])
	fontsize(gca,fs,'points')
	title(nameN(j))
	set(gca,'FontName','Times');


	subplot(2,numNN,(j+numNN))
	%Now use the data in the bivariate distribution to plot curves:
	%Sum up how many CpGs in each class (hypo/hyper/inter)
	Hypobins=find(ygrid<Hypocutoff);
	Hyperbins=find(ygrid>Hypercutoff);
	Interbins=setdiff(1:numel(ygrid),[Hyperbins,Hypobins]);
	
	Hypoarray=sum(Bivariate(:,Hypobins),2);
	Interarray=sum(Bivariate(:,Interbins),2);
	Hyperarray=sum(Bivariate(:,Hyperbins),2);
	%xax=xgrid_P(2:end-1);
	xax=xgridS;
	
	% Normalize the CpG Classifications and plot the curves
	ToCA = Hyperarray+Hypoarray+Interarray;
	%subplot(ax1(j+numel(DataNa)))
	plot(xax,Hypoarray./ToCA,'-o','LineWidth',linsz,'Color','b')
	hold on
	plot(xax,Hyperarray./ToCA,'-o','LineWidth',linsz,'Color','r')
	plot(xax,Interarray./ToCA,'-o','LineWidth',linsz,'Color','k')
	xlabel('Local CpG Density')
	ylabel('Frac. of CpGs')
	xlim([min(xax) max(xax)])
	ylim([0 1])
	ax=gca;
	ax.FontSize=fonttt-4;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];
	set(gca,'FontName','Times');
	title(nameN(j))
	if j==numNN
		legend(['Hypo < ' num2str(Hypocutoff)],['Hyper > ' num2str(Hypercutoff)],'Inter','location','east');
	end
	%Calculating Hypo-Hyper Cross over point
	Hypo=Hypoarray./ToCA;
	Hyper=Hyperarray./ToCA;

	%KB Edit Mar12: Use Finer Gridpoint method to calculate crossover point
	%create a finer grid spacing to calculate the crossover
	finex=linspace(min(xax),max(xax),300);
	%interpolate to a finer grid spacing to calculate the crossover
	Hypofine=interp1(xax,Hypo,finex);
	Hyperfine=interp1(xax,Hyper,finex);
	%find the intersection of the fine curves
	AAA = (Hypofine-Hyperfine);
	[v,inn]=min(abs(AAA));
	title(['ED50 = ' num2str(finex(inn))],'FontSize',(fonttt-4))
	QuanCr(j)=finex(inn);

	IsoCpGinds=find(x==0);
	figure(3000)
	subplot(1,numNN,j)
	histogram(y(IsoCpGinds))
	title(nameN(j))

	MedMeth=[];
	MeanMeth=[];
	Prc25=[];
	Prc75=[];
	for jj=1:numel(densityvals)
    	getx=find(x==densityvals(jj));
    	thesey=y(getx);
    	MedMeth(jj)=median(y(getx));
    	MeanMeth(jj)=mean(y(getx));
    	Prc25(jj)=prctile(thesey,25);
    	Prc75(jj)=prctile(thesey,75);
	end
	title(nameN(j))
	set(gca,'FontName','Times');


	figure(4000)
	subplot(1,numNN,j)
	plot(densityvals,MedMeth,'-o','LineWidth',linsz)
	hold on
	plot(densityvals,MeanMeth,'-o','LineWidth',linsz)
	plot(densityvals,Prc25,'-o','LineWidth',linsz)
	plot(densityvals,Prc75,'-o','LineWidth',linsz)
	xlabel('Local CpG Density')
	ylabel('Methyl. Frac')
	title(nameN(j))
	ax=gca;
	ax.FontSize=fonttt-4;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];
	axis square
	set(gca,'FontName','Times');


	ED50=QuanCr(j);

	%Fitting Hill Function to Ind. CpG Data (Mean)
	choose=MeanMeth;
	normbnds = [min(choose) max(choose)]; 
	Data_1=(choose-normbnds(1))./(normbnds(2)-normbnds(1)); % Normalize the data
	fun = @(p) Fit_Model(p,QuanCr(j),densityvals(:),Data_1(:));

	%Function Minimization metaparameters
	p0=[5]; %initial parameter guess
	mns=[0.005]; %param. minimum 
	mxs=[50]; %param maximum
	options=optimoptions('fmincon','Display','iter');
	nonlcon=[];
	p = fmincon(fun,p0,[],[],[],[],mns,mxs,nonlcon,options);
	H(j)=p; %Save the relevant coefficients
	model_res = ED50.^H(j)./((ED50.^H(j)) + ((densityvals).^H(j)));
	
	%Repeat fitting Hill Function to Ind. CpG Data, using log-transform
	%method of finding the steepest slope.
	%the log transformation method
	newy=log10(Data_1./(1-Data_1));
	newx=log10(densityvals);
	
	%find the halfway point (similar to ED50)
	[vHill,iHill]=min(abs(Data_1-0.5));
	
	%take just the neighboring points in vicinity
	fx=newx(iHill-2:iHill+2);
	fy=newy(iHill-2:iHill+2);
	
	%fit a line to the reduced, log-transformed data
	P = polyfit(fx,fy,1);
	
	%recover the Hill parameters (see equations at top for explanation)
	n=-P(1); 
	K=10^-(P(2)/P(1));
	model_res2=(K^n)./(densityvals.^n+K^n);
	H_med(j)=n; %Save the relevant coefficients for log-transform fitting

	figure(8000)
	set(gca,'FontName','Times');
	subplot(2,numNN,j)
	hold on
	plot(densityvals,Data_1,'ok','LineWidth',linsz*2) % Data
	plot(densityvals,model_res,'-r','LineWidth',linsz) % Model Fit
	if j==numNN
		legend('MeanMeth','Hill Fit')
	end
	if j==1
		ylabel('Norm. Methyl Frac.')
	end
	%xlabel('Local CpG Density')
	title(nameN(j))
	ax=gca;
	ax.FontSize=fonttt-4;
	ax.TickDir='out';
	ax.TickLength=[0.01 0.01];
	axis square
	grid on
	set(gca,'FontName','Times');

	subplot(2,numNN,(j+numNN))
	hold on
	plot(densityvals,Data_1,'ok','LineWidth',linsz*2)
	plot(densityvals,model_res2,'-r','LineWidth',linsz)
	if j==numNN
		legend('MeanMeth','Hill Fit (Log-Tran.)')
	end
	if j==1
		ylabel('Norm. Methyl Frac.')
	end
		xlabel('Local CpG Density')
	%title(nameN(j))
	ax=gca;
	ax.FontSize=fonttt-4;
	axis square
	grid on
	set(gca,'FontName','Times');

end

%% Add in Legends
figure(1000)
subplot(2,numNN,j)
h=colorbar;
ylabel(h,'-log(Probability)')

figure(4000)
legend('MedMeth','MeanMeth','Prc25','Prc75')

%%
figure(8002)
colororder({'#000000','#FF5733'})
yyaxis left
bxx=categorical(nameN);
bxx=reordercats(bxx,nameN);
nil=[nan;nan;nan];
bar(bxx,[H' nil],'grouped')
ylabel('Simple Method')
yticks(0:1:5)

yyaxis right
bar(bxx,[nil H_med'],'grouped')
ylabel('Log-Tran. Method')
yticks(0:5:35)
ax=gca;
ax.FontSize=fonttt;
%axis square
grid on
title('Fitted Hill Coeff')
set(gca,'FontName','Times');

%%
function SSD = Fit_Model(H,ED50,x,Data)
	Hill_fun = ((ED50.^H)./((ED50.^H) + (x.^H)));
	Diff=Data-Hill_fun;
	SSD=sum(Diff.^2); %sum of squared difference, data vs model 
end