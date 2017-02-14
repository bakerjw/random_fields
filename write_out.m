%write_out%%.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%1.04.2010

%This function writes output files for values and locations of the elements
%at each scale

%Calls with: USED_MATRIX(globally)
%Returns:    scale01_output.txt and scale02_output.txt to current directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = write_out()

global USED_MATRIX

%Divide USED_MATRIX into files of seperate scales:
SCALE01 = [];
count01 = 0;
SCALE02 = [];
count02 = 0;
for i = 1:length(USED_MATRIX(:,1))
    if USED_MATRIX(i,7) == 1
        count01 = count01 + 1;
        SCALE01(count01,:) = USED_MATRIX(i,:);
    else
        count02 = count02 + 1;
        SCALE02(count02,:) = USED_MATRIX(i,:);
    end
end
SCALE01(:,2:3) = [];
if count02 > 0
    SCALE02(:,2:3) = [];
end
SCALE01(:,7) = [];
        
%Save seperate scale matrices to seperate output files:
fid = fopen('scale01_output.txt', 'wt');
fprintf(fid, '%6.0f %10.4f %10.4f %10.5f %6.0f %6.0f\n', SCALE01');
fclose(fid);
fid = fopen('scale02_output.txt', 'wt');
fprintf(fid, '%6.0f %10.4f %10.4f %10.5f %6.0f %6.0f %6.0f\n', SCALE02');
fclose(fid);

%end of file