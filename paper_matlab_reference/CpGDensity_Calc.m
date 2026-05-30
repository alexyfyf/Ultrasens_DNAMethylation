%this script outputs an excel file with the following columns:
%'CpGPosition','NumNeighbors','Density'

clear all
filename='chr1_cpgs.csv';
W=1000; %user-defined window size in bp. CpGs within this distance will be considered "neighbors"
WCpG=ceil(W/2); %the number of CpGs to query, on either side of the current CpG
%e.g., if W=50 then we only need to search within 25 CpGs up/down stream
%for potential neighbors
MaxNN=WCpG*2; %the theoretical maximum number of neighbors (if CpGs were directly adjacent)

%load the file
T=readtable(filename);

%extract the CpG positions, using "start" column
Pos=T.start;
%(optional) clear the rest of the data to limit memory usage
clear T

Tot=numel(Pos); %the total number of CpGs in the file
NumKeep=Tot; %the max number of CpGs to keep in the analysis
%To process the whole file, set NumKeep=Tot. setting to a smaller number can
%be useful for debugging/visualization

%determine the starting and ending indices of CpGs to analyze. Two
%considerations: (1) we want to eliminate CpGs too close to beginning or
%end of file--just not sure how to compare their "1-sided" densities with CpGs that
%have neighbors both up/downstream. There are so few of them, eliminating
%them shouldn't affect anything. (2) Keep total number <= NumKeep
stind=WCpG+1; %the starting index to analyze from file
enind=min([NumKeep,Tot-(WCpG)]); %the ending index

allinds=1:NumKeep; %list of indices of the CpGs in the file
allinds=allinds(:);
query=[allinds-WCpG,allinds+WCpG]; %list of query indices [start,end]
keepinds=stind:enind; %indices of CpGs that will be queried
keepPos=Pos(keepinds); %the positions of CpGs kept in the analysis
Num=numel(keepinds); %the number of CpGs kept in the analysis
neighbordists=zeros(Num,2*WCpG+1); %initialize a matrix. Each row is a CpG, %each
%column will contain distances to a potential neighbor. The array size is fixed for code
%efficiency, regardless of actual number of neighbors of each CpG

%loop over the CpGs and determine the neighbors
%although this is a loop, it is somewhat efficient because the array sizes
%remain fixed
for ind=1:Num
    filei=keepinds(ind); %the index of the current CpG in the file
    thiscpg=Pos(filei); %the position of tthis CpG
    query_inds=query(filei,1):query(filei,2); %the file indices to query for neighbors
    query_cpgs=Pos(query_inds); %the CpG positions of query CpGs up/downstream
    dist_cpgs=query_cpgs-thiscpg; %the distances between this CpG and the query set
    neighbor_logical=abs(dist_cpgs)<=W; %find out if they're within distance W
    numneighbors(ind)=sum(neighbor_logical)-1; %the number of neighbors for this CpG 
    %subtract 1 so itself doesn't get counted
    %Create an array: neighbordists. Each row contains the list of distances to neighbors for this CpG.
    %Any potential neighbors that are too far away are given 0 value.
    neighbordists(ind,:)=dist_cpgs.*neighbor_logical;
    %note: this array is not currently being saved or used. But it's
    %created in case it would be useful for further analysis. 
end

%Convert numneighbors to a density using theoretical max number neighbors
Densities=numneighbors/MaxNN;

%NOTE: THIS IS HARD CODED FOR CHR1; WE NEED TO MAITNAIN THE CHROMOSOME
%DESIGNATION FOR BED FILES
ChrDumm=ones(1,length(Densities));
Chr_n='chr1';
Name_n='CpG';
%DataTable=table(keepPos(:),numneighbors(:),Densities(:));
DataTable=table(ChrDumm(:),keepPos(:),[keepPos(:)+1],ChrDumm(:),Densities(:));
%DataTable2=[ChrDumm(:),keepPos(:),[keepPos(:)+1],numneighbors(:),Densities(:)];

%DataTable.Properties.VariableNames={'CpGPosition','NumNeighbors','Density'};
DataTable.Properties.VariableNames={'chr','CpGstart','CpGend','Name','Density'};

% Note that a .bed file requires the following format: chromosome (string),
% start, end, name (optional), score (optional). In our case, the "score"
% would be the density.

DataTable2 = convertvars(DataTable,{'chr','Name'},'string');
DataTable2.chr(:)=deal([Chr_n]);
DataTable2.Name(:)=deal([Name_n]); % Hard-code "names" column to just be CG

outfile=['CpGDensities_W' num2str(W) '.csv'];
%writetable(DataTable,outfile);

outfile2=['CpGDensities_W' num2str(W) '.txt'];
writetable(DataTable2,outfile2,'Delimiter','tab','WriteVariableNames',0);

%Table Head checks
A=DataTable(1:100,:);
AA=DataTable2(1:100,:);