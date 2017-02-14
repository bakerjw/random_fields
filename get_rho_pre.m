%get_rho_pre.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function precalculates correlation versus distance in and between each
%scale.  Precalculation or correlation saves time during the simulation.

%Call with: m, n, ds, d_base, a, b, ROTATE, flag
%Return:    RHO

    %RHO stores correlation versus distance in and between scales.  
    %(:,:,1) is a matrix of possible x distances.
    %(:,:,2) is a matrix of possible y distances.
    %(:,:,3) is a matrix of correlation between elements in scale 2 (fine).
    %(:,:,4) is a matrix of correlation between element in scale 1 and another in scale 2.
    %(:,:,5) is a matrix of correlation between elements in scale 1 (coarse).
      %For a given row and column index (i,j), correlation at x distance (i,j,1)
      %and y distance (i,j,2) is equal to (i,j,[3,4,or 5]).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RHO] = get_rho_pre(m,n,ds,d_base,a,b,ROTATE,flag)

fprintf('\nPrecomputing correlation vs. normalized distance...\n')

%determine h_rel, the maximum normalized distance over which the variogram
%produces non-negligible correlations:
rho_min = 0.0001;  %--> minimum correlation to consider non-negligible

%h_rel is the normalized distance associated with rho_min, and is
%calculated using the inverse of the variogram models.
if flag == 1
    %exponential model
    h_rel = -log(rho_min);
elseif flag == 2
    %Gaussian model
    h_rel = sqrt(-log(rho_min));
elseif flag == 3
    %Gaussian model
    h_rel = 0.596071637983;
end

%determine x_max and y_max from h_rel
ROT_INV = inv(ROTATE);
x_max = max(ROT_INV(1,1)*h_rel*a, ROT_INV(1,2)*h_rel*b);
y_max = max(ROT_INV(2,1)*h_rel*a, ROT_INV(2,2)*h_rel*b);

%set up vectors of possible DX and DY at coarse scale, but limit them to the
%maximums just derived
DX_max = n*d_base;
DY_max = m*d_base;
while(DX_max - d_base > x_max)
    DX_max = DX_max - d_base;
end
while(DY_max - d_base > y_max)
    DY_max = DY_max - d_base;
end
DX1 = -DX_max:d_base:DX_max;
DY1 = -DY_max:d_base:DY_max;

%Get denominator for correlation calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x&y modification vectors (Xmod and Ymod) will determine fine scale (x,y) 
%coordinates from coarse scale starting point
DS = 1*ds;
increment12 = d_base/DS;
startx = increment12*(1-DS)/2;
starty = increment12*(1-DS)/2;
Xrow = startx:increment12:(startx + increment12*(DS-1));
Yrow = starty:increment12:(starty + increment12*(DS-1));
[Xgrid,Ygrid] = meshgrid(Xrow,Yrow);
Xmod = Xgrid(:);
Ymod = Ygrid(:);

%set reference point.  Use (x,y)=(0,0).
x1 = 0;
y1 = 0;

%Get smallxy1, a vector of x,y coordinates of the fine-scale elements within
%a coarse-scale element, then determine the possible x,y distances and 
%associated correlations between them.
bigxy1 = [x1,y1];
smallxy1 = [bigxy1(1)+Xmod,bigxy1(2)+Ymod];
[XA,XB] = meshgrid(smallxy1(:,1),smallxy1(:,1));
D1 = XA - XB;
D1 = D1(:);
[YA,YB] = meshgrid(smallxy1(:,2),smallxy1(:,2));
D2 = YA - YB;
D2 = D2(:);

rho = variogram(D1,D2,a,b,ROTATE,flag);
RHO_denominator = sum(rho);   %to be used for correlation calculations involving
  %coarse-scale elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set up x and y modification vectors at an even finer resolution than the
%fine scale.  This step is for the eventual creation of an array of x,y
%distance pairs for which to precalculate correlation, and in order to
%better see trends in the plots the resolution is increased beyond the fine
%scale.  This comes at a minor computational cost.
DS = 1*ds;
increment12 = d_base/DS;
startx = increment12*(1-DS)/2;
starty = increment12*(1-DS)/2;
Xrow = startx-increment12/2:increment12/2:(startx + increment12*(DS-1)+increment12/2);
Yrow = starty-increment12/2:increment12/2:(starty + increment12*(DS-1)+increment12/2);
[Xgrid,Ygrid] = meshgrid(Xrow,Yrow);
DX2 = Xgrid(:); %similar to Xmod above, but for finer scale
DY2 = Ygrid(:); %similar to Ymod above, but for finer scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get vectors of all possible dx's and dy's
[XA,XB] = meshgrid(DX1,DX2);
Xdist = XA - XB;
Xdist = Xdist(:);
Xdist = sortrows(Xdist);
%reduce size to include only one repetition of each distance
g = length(Xdist);
i = 2;
if g >= i
    while g >= i
        if abs(Xdist(i) - Xdist(i-1)) < 0.0001
            Xdist(i) = [];
            g = length(Xdist);
        else
            i = i + 1;
        end
    end
end
[YA,YB] = meshgrid(DY1,DY2);
Ydist = YA - YB;
Ydist = Ydist(:);
Ydist = sortrows(Ydist);
%reduce size to include only one repetition of each distance
g = length(Ydist);
i = 2;
if g >= i
    while g >= i
        if abs(Ydist(i) - Ydist(i-1)) < 0.0001
            Ydist(i) = [];
            g = length(Ydist);
        else
            i = i + 1;
        end
    end
end

%reduce dx's and dy's to within maximum limits
i = length(Xdist);
while(Xdist(i) > x_max)
    Xdist(i) = [];
    Xdist(1) = [];
    i = length(Xdist);
end
y_cut = length(Ydist) - round(length(Ydist)/2);
Ydist(1:y_cut) = [];
i = length(Ydist);
while(Ydist(i) > y_max)
    Ydist(i) = [];
    i = length(Ydist);
end

%get correlations using variogram.m
for g = 1:length(Xdist)
    for i = 1:length(Ydist)
        xplace = g;
        yplace = i;
        RHO(yplace,xplace,1) = Xdist(g);
        RHO(yplace,xplace,2) = Ydist(i);
        
        bigxy2 = [Xdist(g),Ydist(i)];
        smallxy2 = [bigxy2(1)+Xmod,bigxy2(2)+Ymod];
                
        %Calculate all rho_fine_fine in one step
        D1 = bigxy2(1) - x1;
        D2 = bigxy2(2) - y1;
        RHO(yplace,xplace,3) = variogram(D1,D2,a,b,ROTATE,flag);
        
        %get rho_avg_avg
        [XA,XB] = meshgrid(smallxy1(:,1),smallxy2(:,1));
        D1 = XB - XA;
        D1 = D1(:);
        [YA,YB] = meshgrid(smallxy1(:,2),smallxy2(:,2));
        D2 = YB - YA;
        D2 = D2(:);
        rho = variogram(D1,D2,a,b,ROTATE,flag);
        RHO_numerator = sum(rho);
        RHO(yplace,xplace,5) = RHO_numerator/RHO_denominator;
        
        %get rho_fine_avg
        ONES = ones(DS*DS,1);
        D1 = smallxy2(:,1) - ONES*x1;
        D2 = smallxy2(:,2) - ONES*y1;
        rho = variogram(D1,D2,a,b,ROTATE,flag);
        RHO_numerator = sum(rho);
        RHO(yplace,xplace,4) = RHO_numerator/sqrt(RHO_denominator);
    end
end

%PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotflag = 1; % 1 displays plots, any other value will not.
if plotflag==1
    Lm = length(RHO(:,1,1));
    Ln = length(RHO(1,:,1));
    zerom = Lm - round(Lm/2) + 1;
    zeron = Ln - round(Ln/2) + 1;
    figure(10)
    hold on
    plot(RHO(:,zeron,2),RHO(:,zeron,5),':k', 'linewidth',1.5)
    plot(RHO(:,zeron,2),RHO(:,zeron,4),'--b','linewidth',1.5)
    plot(RHO(:,zeron,2),RHO(:,zeron,3),'-r', 'linewidth',1.5)
    set(gca, 'Fontsize', 12)
    legend('\rho_{Z1,Z1}','\rho_{Z1,Z2}','\rho_{Z2,Z2}');
    xlabel('Normalized separation distance, h', 'Fontsize', 14);
    ylabel('\rho', 'Fontsize', 14);
    xlim([0 4.5])
    legend boxoff
    hold off
    box off
end

%end of file