%get_var_avg.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function calculates the variance for coarse-scale elements

%Call with: ds, d_base, sigma, a ,b ,ROTATE, flag
%Return:    rho

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [var_avg] = get_var_avg(ds,d_base,sigma,a,b,ROTATE,flag)

scale_factor = 2;       %always want SF to be 2 for this function
bigxy = [0.5,0.5];      %dummy value for coarse-scale location at which to calculate variance

for i = 1:ds*ds                                                            %this loop attaches i,j,x,y to indices
    index_i = ceil(i/ds);
    if i-(index_i-1)*ds==0			
        index_j = ds;
    else
        index_j = i-(index_i-1)*ds;
    end
    index_i = ds - index_i + 1;                                            %reverses index numbering scheme. Turn off to return to default
    index_x = bigxy(1,1) - d_base*(1/(2*ds^(scale_factor-1)))*(1+ds-2*index_j);
    index_y = bigxy(1,2) + d_base*(1/(2*ds^(scale_factor-1)))*(1+ds-2*index_i);
    smallxy(i,1) = index_x;
    smallxy(i,2) = index_y;
end

rho_var_sum = 0;

for i = 1:ds*ds
    for j = 1:ds*ds
        d1 = smallxy(i,1)-smallxy(j,1);
        d2 = smallxy(i,2)-smallxy(j,2);
        rho = variogram(d1,d2,a,b,ROTATE,flag);
        rho_var_sum = rho_var_sum + rho*sigma^2;
    end
end

var_avg = rho_var_sum / ds^4;

%end of file