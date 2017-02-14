%get_indices.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function determines spatial x and y coordinates as well as local i,j 
%row,column coordinates.

%Call with: m, n, ds, d_base, scale_factor, x_big, y_big
%Return:    indices

    %indices(1,1) is index (local) of next point.  Each element in the
      %field is assigned a number, or "index".
    %indices(2,1) is i location (local) in grid of next point.
    %indices(3,1) is j location (local) in grid of next point.
    %indices(4,1) is x location (global) in grid of next point.
    %indices(5,1) is y location (global) in grid of next point.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indices = get_indices(m, n, ds, d_base, scale_factor, x_big, y_big)

if scale_factor == 1
    %create indices for coarse-scale field
    [temp,indices] = sort(unifrnd(0,1,[1,m*n]));                           %creates randomized vector of indices
    for i = 1:m*n                                                          %this loop attaches i,j,x,y to indices
        index_i = ceil(indices(1,i)/n);
        if indices(1,i)-(index_i-1)*n==0			
            index_j = n;
        else
            index_j = indices(1,i)-(index_i-1)*n;
        end
        index_i = m - index_i + 1;                                         %reverses index numbering scheme. Turn off to return to default
        index_x = index_j*d_base - d_base/2;
        index_y = d_base*m-(index_i*d_base - d_base/2);
        indices(2,i) = index_i;
        indices(3,i) = index_j;
        indices(4,i) = index_x;
        indices(5,i) = index_y;
    end
else
    %create indicies of fine-scale elements within a coarse-scale element
    [temp,indices] = sort(unifrnd(0,1,[1,ds*ds]));                         %creates randomized vector of indices
    for i = 1:ds*ds                                                        %this loop attaches i,j,x,y to indices
        index_i = ceil(indices(1,i)/ds);
        if indices(1,i)-(index_i-1)*ds==0			
            index_j = ds;
        else
            index_j = indices(1,i)-(index_i-1)*ds;
        end
        index_i = ds - index_i + 1;                                        %reverses index numbering scheme. Turn off to return to default
        index_x = x_big - d_base*(1/(2*ds^(scale_factor-1)))*(1+ds-2*index_j);
        index_y = y_big + d_base*(1/(2*ds^(scale_factor-1)))*(1+ds-2*index_i);
        indices(2,i) = index_i;
        indices(3,i) = index_j;
        indices(4,i) = index_x;
        indices(5,i) = index_y;
    end
end

%end of file