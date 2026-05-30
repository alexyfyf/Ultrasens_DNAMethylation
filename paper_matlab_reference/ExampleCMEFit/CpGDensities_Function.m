function Densities = CpGDensities_Function(CpGPositions,W)

%inputs
%CpGPositions - vector of positions in bp
%W window size, it counts neighbors w/in +/- W


NCpG=length(CpGPositions);

WCpG=ceil(W/2); %the number of CpGs to query, on either side of the current CpG
%e.g., if W=50 then we only need to search within 25 CpGs up/down stream
%for potential neighbors
MaxNN=WCpG*2; %the theoretical maximum number of neighbors (if CpGs were directly adjacent)

%loop over the CpGs and determine the neighbors
%although this is a loop, it is somewhat efficient because the array sizes
%remain fixed
for ind=1:NCpG
    thiscpg=CpGPositions(ind); %the position of tthis CpG
    query_inds=1:NCpG;
    query_cpgs=CpGPositions(query_inds);
   
    dist_cpgs=query_cpgs-thiscpg; %the distances between this CpG and the query set
    neighbor_logical=abs(dist_cpgs)<=W; %find out if they're within distance W
    numneighbors(ind)=sum(neighbor_logical)-1; %the number of neighbors for this CpG 
    %subtract 1 so itself doesn't get counted
end

%Convert numneighbors to a density using theoretical max number neighbors
Densities=numneighbors/MaxNN;

end