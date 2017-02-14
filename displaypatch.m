%displaypatch.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function plots the simulated sample.

%Call with: d_base, ds
%Return:    visual display of the sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = displaypatch(d_base, ds, coarseFlag)

global USED_MATRIX
global HISTORY

if coarseFlag == 1                   %Plot refined field
    WORKING = USED_MATRIX;
else                                 %Plot coarse field
    indices = HISTORY(:,7) == 1;
    WORKING = HISTORY(indices,:);
end

delete_row = 'temp';
while isempty(delete_row) ~= 1
    delete_row = find(WORKING(:,8) == 1,1,'first');
    WORKING(delete_row,:) = [];
end

Z_max = max(max(WORKING(:,6)));
Z_min = min(min(WORKING(:,6)));
X = [];
Y = [];

for i = 1:length(WORKING(:,1))
    scale_factor = WORKING(i,7);
    if scale_factor == 1
        X = [X [WORKING(i,4)-d_base/2;WORKING(i,4)-d_base/2;WORKING(i,4)+d_base/2;WORKING(i,4)+d_base/2]];
        Y = [Y [WORKING(i,5)-d_base/2;WORKING(i,5)+d_base/2;WORKING(i,5)+d_base/2;WORKING(i,5)-d_base/2]];
    else
        X = [X [WORKING(i,4)-d_base/(2*ds^(scale_factor-1));WORKING(i,4)-d_base/(2*ds^(scale_factor-1));WORKING(i,4)+d_base/(2*ds^(scale_factor-1));WORKING(i,4)+d_base/(2*ds^(scale_factor-1))]];
        Y = [Y [WORKING(i,5)-d_base/(2*ds^(scale_factor-1));WORKING(i,5)+d_base/(2*ds^(scale_factor-1));WORKING(i,5)+d_base/(2*ds^(scale_factor-1));WORKING(i,5)-d_base/(2*ds^(scale_factor-1))]];
    end
end


p = patch(X,Y,ones(size(X)) ,'b');
clear cdata
set(gca,'CLim',[Z_min Z_max])
set(p,'FaceColor','flat','FaceVertexCData',WORKING(:,6),'CDataMapping','scaled')
colormap jet
colorbar
