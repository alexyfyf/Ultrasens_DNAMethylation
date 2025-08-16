function [n,K,PlotModel,err] = FitHill2Param(x,y,Direction)
%Direction says whether the function is decreasing (0, e.g., Mean) or
%increasing (1, e.g., Hypo curve)
basex=x;
%remove the stuff at high density (noisy)
 % getinds=x<0.4;
 % x=x(getinds);
 % y=y(getinds);

%interpolate to get more fine curve
finex=linspace(min(x),max(x),300);
finey=interp1(x,y,finex);

x=finex;
y=finey;

x=x(:);
%correct offset from 0 and normalize
%Miny=min(y);
%y=y-Miny;

[Maxy,maxi]=max(y);
y=y/Maxy;
%remove the part at the beginning, if it dips
%y=y(maxi:end);
%x=x(maxi:end);

Data=y(:);
%Fitting functioni
fun =@(p)Fit_Model(p,Data,x,Direction);
p0=[1.1,1.1]; %initial parameter guess [K,n]
mns=[0,0]; %parameter range: minimum
mxs=[20,20]; %maximum

options=optimoptions('fmincon','Display','none');%,'Display','iter')
nonlcon=[];
%Perform the fitting:
p = fmincon(fun,p0,[],[],[],[],mns,mxs,nonlcon,options);

K=p(1);
n=p(2);

if Direction
    PlotModel=Maxy*(basex.^n)./(K^n+basex.^n); %model to be fit
else
    PlotModel=Maxy*(K^n)./(K^n+basex.^n); %model to be fit
end

err=Fit_Model(p,Data,x,Direction);

    function SSD = Fit_Model(p,Data,x,Direction)
        K=p(1);
        n=p(2);

        if Direction
            Model=(x.^n)./(K^n+x.^n); %model to be fit
        else
            Model=(K^n)./(K^n+x.^n); %model to be fit
        end

        Diff=Data-Model;
        %remove the upper part of curve where data is noisy
        SSD=sum(Diff.^2); %sum of squared difference, data vs model
    end
end

