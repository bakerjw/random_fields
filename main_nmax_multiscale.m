%MAIN FILE: MULTI-SCALE RANDOM FIELD SIMULATION SCRIPT

%Andrew Seifried
%Stanford University
%12.27.2009
%Updated 12.22.2010

%This script generates a two-dimensional multi-scale simulation of a random field in
%standard normal space and transforms those values into a specified distribution.  
%The framework is currently set up to use two scales.  A coarse scale field of
%random variables is initially simulated, then some or all of those elements may be 
%subdivided to include one finer scale.  

%The user may define in advance which coarse-scale elements to subdivide,
%or can specify a critical value of either the difference in value of neighboring
%simulated elements or the variance of its neighbors.
%The first option is done by specifying the (x,y) coordinates of the center of
%the coarse-scale element to be subdivided in the text file 'subdivide.txt',
%which needs to be included in the current directory.  The other options
%are completed by modifying input values below.

%The user must also specify the target distribution and incorporate it into
%'transform.m'.

%A text file for each scale is saved in the current directory:
%scale01_output.txt and scale02_output.txt.

%Variables and functions are generally described as they are introduced.
%See the end of this script for more detailed descriptions of the more complex 
%variables used and all of the m files created and used here.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear;clc;

%GLOBAL DECLARATIONS (list of all global variables in the program)
global indices           %stores x,y (spatial location) and i,j (matrix row,column)
                         %coordinates of each element.
global USED_MATRIX       %stores information only regarding current refined state 
                         %(as elements are refined, coarser-scale simulated values 
                         %are removed from the array).
global HISTORY           %stores information of all previously simulated elements.
global RHO               %stores calculated values of correlation versus distance.


indices=[];USED_MATRIX=[];n_max=[];mu=[];sigma=[];sigma_avg=[];HISTORY=[];RHO=[];

%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SOFT INPUT VALUES (these are specified by the user)
m = 	10;         %Number of elements in each column of coarse grid.
n = 	5;          %Number of elements in each row of coarse grid.
ds =    4;          %ds = "division size".  Coarse elements are subdivided into ds x ds finer elements. ds MUST be even!
d_base =1;          %Units of distance comprising 1 side of coarsest element.
theta = 45;         %Angle of principal axis (+ deg. CCW from horizontal).
a =     10;          %Dictates how quickly correlation decays with distance along the principal axis.
b =     1;          %Dictates how quickly correlation decays with distance perpendicular to the principal axis.
n_max = 125;        %Maximum number of previously generated elements to correlate to.
variogramFlag = 1;  %Determines which variogram to use
                        %variogramFlag == 1 --> exponential model
                        %variogramFlag == 2 --> Gaussian model
                        %variogramFlag == 3 --> spherical model
cutoff = 1.0;       %Cutoff value used as critical value in chosen subdivision method.
cutoffFlag = 1;     %Determines which subdivision method to use
                        %cutoffFlag == 0 --> subdivision based only on the text file "subdivide.txt"
                            % **the text file is also used for the remaining options**
                        %cutoffFlag == 1 --> subdivision based on the maximum difference with any neighbor
                        %cutoffFlag == 2 --> subdivision based on the variance of neighboring values
                        %cutoffFlag == 3 --> subdivision based on the variance of neighboring values, includes
                            %more neighboring elements than previous variance-based option.
                        
%HARD INPUT VALUES (do not change these)
mu =	0;          %Mean of standard normal RV at finest scale (scale 2 in two-scale example).
sigma =	1;          %Standard deviation of standard normal RV at finest scale.
scale_factor = 1;   %Initial scale factor is always 1 (simulation begins isn the coarse scale)

%PRECALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Precalculate Rotation Matrix
rad = theta*pi/180;
ROTATE = [cos(rad), sin(rad); -sin(rad), cos(rad)];

%Precalculate standard deviation of the coarse scale:
fprintf('\nPrecomputing standard deviation of coarse scale...\n')
var_avg = get_var_avg(ds,d_base,sigma,a,b,ROTATE,variogramFlag);
sigma_avg = sqrt(var_avg);

%Precalculate possible combinations of rho in and between each scale
RHO = get_rho_pre(m,n,ds,d_base,a,b,ROTATE,variogramFlag);

%Set up first set of indices (at coarse scale)
indices = get_indices(m, n, ds, d_base, scale_factor);                     %randomized indices of points to generate, determines x and y coordinates

SIGS = [];          %records conditional sigmas as they are calculated
sig_count1 = 1;     %i position of SIGS
sig_count2 = 1;     %j position of SIGS

%SIMULATE COARSE-SCALE FIELD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fprintf('\nSimulating coarse field...\n')

Z=zeros(m,n);  % --> 'Z' denotes standard normal RV, 'U' denotes RV after transformation to target distribution

i1 = indices(2,1);  %i of first point
j1 = indices(3,1);  %j of first point

Z(i1,j1) = normrnd(mu,sigma_avg);                                          %simulate first element from standard normal distribution
USED_MATRIX = [indices(:,1)', Z(i1,j1), scale_factor,0,0];                 %moves index just used from indices to used
indices(:,1) = [];                                                         %removes index just used from "indices"

index_big = 0;                                                             %no associated up-scaled element at coarse scale
while isempty(indices)~=1
    i1 = indices(2,1);  %i of next point
    j1 = indices(3,1);  %j of next point

    [mu_prime,sigma_prime] = get_stats(scale_factor,a,b,ROTATE,index_big,sigma_avg,n_max,mu,sigma);
    SIGS(sig_count1,sig_count2) = sigma_prime;
    sig_count1 = sig_count1 + 1;
    
    Z(i1,j1) = normrnd(mu_prime,sigma_prime);                              %simulate the next value
    used_input = [indices(:,1)',Z(i1,j1),scale_factor,0,0];
    USED_MATRIX = [USED_MATRIX; used_input];
    indices(:,1) = [];
end
HISTORY = USED_MATRIX;
sig_count2 = sig_count2 + 1;
sig_count1 = 1;

%transform(1,ds,d_base,a,b,ROTATE,variogramFlag,sigma_avg);                %transform from z to u to analyze for subdivisions
subdivide(d_base, cutoff, cutoffFlag);                                    %flags coarse-scale elements for subdivision
%transform(2,ds,d_base,a,b,ROTATE,variogramFlag,sigma_avg);                %transform from u to z to continue simulation


%SUBDIVISION TO FINE-SCALE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nRefining simulated field...\n')
    
scale_factor = scale_factor + 1;

subdiv_rows = find(USED_MATRIX(:,8) == 1, 1, 'first');

while isempty(subdiv_rows) ~= 1
    x_big = USED_MATRIX(subdiv_rows(1),4);
    y_big = USED_MATRIX(subdiv_rows(1),5);
    index_big = USED_MATRIX(subdiv_rows(1),1);

    Z = zeros(ds,ds);

    indices = get_indices(m, n, ds, d_base, scale_factor, x_big, y_big);

    while isempty(indices) ~= 1
        i1 = indices(2,1);  %i of next point
        j1 = indices(3,1);  %j of next point

        [mu_prime,sigma_prime] = get_stats(scale_factor,a,b,ROTATE,index_big,sigma_avg,n_max,mu,sigma);
        SIGS(sig_count1,sig_count2) = sigma_prime;
        sig_count1 = sig_count1 + 1;

        Z(i1,j1) = normrnd(mu_prime,sigma_prime);                          %simulates next fine-scale element
        used_input = [indices(:,1)',Z(i1,j1),scale_factor,0,index_big];
        USED_MATRIX = [USED_MATRIX; used_input];
        HISTORY = [HISTORY; used_input];
        indices(:,1) = [];
    end
    subdiv_rows(1) = [];
    g = find(((USED_MATRIX(:,4) == x_big) & (USED_MATRIX(:,5) == y_big)) & USED_MATRIX(:,7)==scale_factor-1);
    USED_MATRIX(g,:) = [];                                                 %remove larger-scale simulated value of refined element
    subdiv_rows = find(USED_MATRIX(:,8) == 1, 1, 'first');
    sig_count2 = sig_count2 + 1;
    sig_count1 = 1;
end

%END SIMULATION


%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nEnd of simulation.  Transforming to target distribution...\n')

%Uncomment lines below to see histograms of simulated values, if desired
figure
hist(USED_MATRIX(:,6))
title('Histogram of standard normal realizations')

transform(1,ds,d_base,a,b,ROTATE,variogramFlag,sigma_avg);                          %transform from z to u
figure
hist(USED_MATRIX(:,6))
title('Histogram of transformed realizations')

%Write output files
write_out();

%Display the generated field (both the coarse and refined fields)
figure
displaypatch(d_base, ds, 0);
figure
displaypatch(d_base, ds, 1);

fprintf('\nDone.\n')

%Variable definitions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %'n_max' is the maximum number of previously simulated elements to 
      %correlate to.
      
    %'mu' is the mean of an element at any scale.
    
    %'sigma' is the standard deviation of a fine-scale element.
    
    %'sigma_avg' is the standard deviation of a coarse-scale element.
    
    %'indices' keeps track of the location of each element locally.
    %indices(1,1) is index (local) of next point.  Each element in the
      %field is assigned a number, or "index".
    %indices(2,1) is i location (local) in grid of next point.
    %indices(3,1) is j location (local) in grid of next point.
    %indices(4,1) is x location (global) in grid of next point.
    %indices(5,1) is y location (global) in grid of next point.
    
    %As realizations are created for each point, they are stored in 
      %'USED_MATRIX'
    %USED_MATRIX(:,1) is index of point (local).
    %USED_MATRIX(:,2) is i location (local).
    %USED_MATRIX(:,3) is j location (local).
    %USED_MATRIX(:,4) is x location (global).
    %USED_MATRIX(:,5) is y location (global).
    %USED_MATRIX(:,6) is realization of RV at this point.
    %USED_MATRIX(:,7) is scale factor.
    %USED_MATRIX(:,8) is flag to determine subdivision.
    %USED_MATRIX(:,9) is associated element index at previous scale.
    
    %'HISTORY' is structured the same as USED_MATRIX, but after coarse scale
    %elements (scale 1) are refined to the scale 2 elements of which they
    %are comprised, the scale 1 element is deleted.  USED_MATRIX is used
    %for correlation of new elements to previously generated elements and
    %deleting scale 1 avoids "double-counting" elements representing the
    %same area.
    
    %'RHO' stores correlation versus distance in and between scales.  
    %(:,:,1) is a matrix of possible x distances.
    %(:,:,2) is a matrix of possible y distances.
    %(:,:,3) is a matrix of correlation between elements in scale 2 (fine).
    %(:,:,4) is a matrix of correlation between element in scale 1 and another in scale 2.
    %(:,:,5) is a matrix of correlation between elements in scale 1 (coarse).
      %For a given row and column index (i,j), correlation at x distance (i,j,1)
      %and y distance (i,j,2) is equal to (i,j,[3,4,or 5]).

    %'breakpoints' is vertical array of (x,y) locations of points to be
      %subdivided.
    
    %'Z' is a holding matrix for the current grid (m x m or ds x ds) being
      %generated.


%Functions called%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %get_indices.m   returns local i (row) and j (column) and global x and
                     %y (spatial location) values
    %get_rho_pre.m   precalculates correlation versus dustance to save
                     %computational expense.
    %get_stats.m     determines conditional mu and sigma based on
                     %previously simulated elements.
    %variogram.m     determines correlation between points.  Variogram can
                     %handle individual points or arrays of coordinates.
    %get_var_avg.m   determines variance of a coarse-scale element
    %transform.m     transforms standard normal field to target dist'n and
                     %rotates the principal axis
    %subdivide.m     reads text file to flag elements to be subdivided.
                     %Also contains framework to check for either large
                     %differences in simulated values of adjacent points or
                     %large variances of neighboring coarse-scale elements
                     %to flag them for subdivision.
    %displaypatch.m  returns graphical representation of final simulated
                     %multi-scale sample.
    
%end of file