%subdivivide.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%Adaptive refinement methods added by:
%Jack Baker, Stanford University,
%6/2010

%This function flags coarse-scale elements to be subdivided down to the
%finer scale.

%The user must specify the (x,y) coordinates of the center of the
%coarse-scale elements which are to be subdivided in 'subdivide.txt' which
%should be in the the current directory.

%The framework for comparing the values of adjacent simulated elements and
%comparing them to a cutoff value, therby flagging the elements on either
%side of a large gradient for subdivision to achieve more resolution there,
%is included but has not been tested or implemented yet.

%Calls with: USED_MATRIX(globally), d_base, cutoff
%Returns: USED_MATRIX(globally) with updated flag to subdivide

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = subdivide(d_base, cutoff, checkFlag)

global USED_MATRIX

%FIRST, FLAG SPECIFIED POINTS IN THE TEXT FILE:

file_in = fopen('subdivide.txt', 'r');
temp = fgetl(file_in); % read and discard the header line
breakpoints = fscanf(file_in, '%f %f', [2,inf]); % read to end of file
breakpoints = breakpoints';
fclose(file_in);
    
while isempty(breakpoints)~=1
    x = breakpoints(1,1);
    y = breakpoints(1,2);

    flag_row = find((round(USED_MATRIX(:,4)*10e4)/10e4 == x) & (round(USED_MATRIX(:,5)*10e4)/10e4 == y),1);
    USED_MATRIX(flag_row,8) = 1; %sets the subdivide flag to 1

    breakpoints(1,:) = [];
end

%SECOND, CHECK FOR JUMPS IN VALUE OF NEIGHBORING ELEMENTS.

if checkFlag == 1  %DIFFERENCE IN ADJACENT ELEMENT VALUES
    for i = 1:length(USED_MATRIX(:,1))
        CHECK = [0,1;-1,0;1,0;0,-1]; %Checks cells adjacent to current cell
        x = USED_MATRIX(i,4);
        y = USED_MATRIX(i,5);
        POR = USED_MATRIX(i,6);
        while ~isempty(CHECK) 
            temp_row = find((USED_MATRIX(:,4) == x + d_base*CHECK(1,1)) & (USED_MATRIX(:,5) == y + d_base*CHECK(1,2)));
            if ~isempty(temp_row)
                if abs(USED_MATRIX(temp_row,6) - POR) > cutoff
                    %[i POR  temp_row USED_MATRIX(temp_row,6) abs(USED_MATRIX(temp_row,6) - POR)] %output comparison numbers
                    USED_MATRIX(temp_row,8) = 1; %sets flag for subdivision
                    USED_MATRIX(i,8) = 1;
                end
            end
            CHECK(1,:) = [];
        end
    end
    
elseif checkFlag == 2  %CHECK VARIANCE OF NEIGHBORING ELEMENTS
    for i = 1:length(USED_MATRIX(:,1))
        simulated(USED_MATRIX(i,2), USED_MATRIX(i,3)) = USED_MATRIX(i,6);
    end
    [m n] = size(simulated);
    
    for k = 1:length(USED_MATRIX(:,1))
        i_vals = [max(USED_MATRIX(k,2)-1,1):min(USED_MATRIX(k,2)+1,m)];
        j_vals = [max(USED_MATRIX(k,3)-1,1):min(USED_MATRIX(k,3)+1,n)];
        neightbors = (simulated(i_vals,j_vals));
        %var(neightbors(:))
    
        if var(neightbors(:)) > cutoff
            USED_MATRIX(k,8) = 1; %sets flag for subdivision
        end
    end
    
elseif checkFlag == 3  %CHECK VARIANCE OF NEIGHBORING ELEMENTS - Include more elements
    for i = 1:length(USED_MATRIX(:,1))
        simulated(USED_MATRIX(i,2), USED_MATRIX(i,3)) = USED_MATRIX(i,6);
    end
    [m n] = size(simulated);
    
    for k = 1:length(USED_MATRIX(:,1))
        i_vals = [max(USED_MATRIX(k,2)-2,1):min(USED_MATRIX(k,2)+2,m)];
        j_vals = [max(USED_MATRIX(k,3)-2,1):min(USED_MATRIX(k,3)+2,n)];
        neightbors = (simulated(i_vals,j_vals));
        %var(neightbors(:))
    
        if var(neightbors(:)) > cutoff
            USED_MATRIX(k,8) = 1; %sets flag for subdivision
        end
    end
end

num_to_refine = sum(USED_MATRIX(:,8))
%end of file