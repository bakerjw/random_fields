%variogram.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function uses a specified variogram model to determine correlation
%between elements (either individually or in arrays).

%Call with: d1, d2, a, b, ROTATE, flag
%Return:    rho

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho] = variogram(d1,d2,a,b,ROTATE,flag)

h1 = ROTATE(1,1)*d1 + ROTATE(1,2)*d2;
h2 = ROTATE(2,1)*d1 + ROTATE(2,2)*d2;
    
h = ((h1./a).^2 + (h2./b).^2).^(1/2); %variables 'a' and 'b' determine how quickly correlation dies with distance.

if flag == 1
%Exponential model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = 1-exp(-h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif flag == 2
%Gaussian model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = 1-exp(-h.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

elseif flag == 3
%Spherical model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = (3/2).*h - (1/2).*h.^3;
    r = find(h > 1);
    gamma(r) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

rho = 1 - gamma;

%end of file