%get_stats.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function determines the conditional mean and standard deviation based
%on previously simulated values.

%Call with: scale_factor, a, b, ROTATE, index_big
%Return:    mu_prime, sigma_prime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu_prime, sigma_prime] = get_stats(scale_factor, a, b, ROTATE,...
    index_big, sigma_avg, n_max, mu,sigma)

global indices
global USED_MATRIX
global RHO

%calculate covariance matrix when FEWER than nmax elements have been
  %previously simulate:
l = length(USED_MATRIX(:,1));
if l <= n_max
        
    SIGMA12 = zeros(1,l);
    SIGMA22 = zeros(l,l);
           
    for i = 1:l
        
        %Get next entry in SIGMA12, the covariance between the new element
        %and previously simulated elements.
        X1 = USED_MATRIX(i,4);      %x location of first element
        Y1 = USED_MATRIX(i,5);      %y location of first element
        SF1 = USED_MATRIX(i,7);     %scale factor of first element
        X2 = indices(4,1);          %x location of second element
        Y2 = indices(5,1);          %y location of second element
        SF2 = scale_factor;         %scale factor of second element
        
        dX = X2 - X1;               %x distance
        dY = Y2 - Y1;               %y distance
        if dY~=0
            dX = sign(dX./dY) * abs(dX);
        end
        dY = abs(dY);
        
        scale_index = (SF1 + SF2) - 1;
        
        %interpolate correlation from previously calculated RHO matrix.
        if SF2 == 1
            rho11 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,5), dX, dY,'*linear');
            rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
            rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
            rho_temp = [rho11,rho12,rho22];
            [r,c] = find(isnan(rho_temp));
            rho_temp(r,c) = 0;
            rho = rho_temp(:,scale_index);   
        else
            rho11 = 0;
            rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
            rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
            rho_temp = [rho11,rho12,rho22];
            [r,c] = find(isnan(rho_temp));
            rho_temp(r,c) = 0;
            rho = rho_temp(:,scale_index);
        end
        
        sigma_temp = [sigma_avg, sigma];
        sigma_n = sigma_temp(:,scale_factor);
        sigma_1 = sigma_temp(:,SF1);
        SIGMA12(1,i) = rho*sigma_n*sigma_1;
        
        %Get next row of SIGMA22, the covariances between previously
        %simulated elements.
        X1 = X1 * ones(l,1);
        Y1 = Y1 * ones(l,1);
        SF1 = SF1 * ones(l,1);
        X2 = USED_MATRIX(:,4);
        Y2 = USED_MATRIX(:,5);
        SF2 = USED_MATRIX(:,7);
    
        dX = X2 - X1;
        dY = Y2 - Y1;
        for p = 1:l
            if dY(p)~=0
                dX(p) = sign(dX(p)./dY(p)) * abs(dX(p));
            end
            dY(p) = abs(dY(p));
        end

        scale_index = (SF1 + SF2) - 1;                   %vector of SF's
        
        %interpolate correlation from previously calculated RHO matrix.
        rho11 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,5), dX, dY,'*linear');
        rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
        rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
        rho_temp = [rho11,rho12,rho22];
        [r,c] = find(isnan(rho_temp));
        rho_temp(r,c) = 0;
        rho = diag(rho_temp(:,scale_index));
        
        sigma_temp = ones(l,1) * [sigma_avg, sigma];
        sigma_1 = ones(l,1) * sigma_1;
        sigma_2 = diag(sigma_temp(:,SF2));
        SIGMA22(i,1:l) = [rho.*sigma_1.*sigma_2]';
    end
else
    
%calculate covariance matrix when MORE than nmax elements have been
  %previously simulate:
    SIGMA12 = zeros(1, n_max);
    SIGMA22 = zeros(n_max, n_max);
    
    %'X' is a variable similar to 'USED_MATRIX', but it gets sorted by
    %distance and then truncated after n_max elements are included.  This
    %way the n_max most relevant elements are considered.
    d1 = USED_MATRIX(:,4)-ones(l,1)*indices(4,1);
    d2 = USED_MATRIX(:,5)-ones(l,1)*indices(5,1);
    h1 = ROTATE(1,1)*d1 + ROTATE(1,2)*d2;
    h2 = ROTATE(2,1)*d1 + ROTATE(2,2)*d2;
    X(:,1) = ((h1./a).^2 + (h2./b).^2).^(1/2);
    X(:,2) = USED_MATRIX(:,4);
    X(:,3) = USED_MATRIX(:,5);
    X(:,4) = USED_MATRIX(:,6);
    X(:,5) = USED_MATRIX(:,7);
    X(:,6) = USED_MATRIX(:,1);
    X(:,7) = USED_MATRIX(:,9);
    
    X = sortrows(X);
    
    %Also include all elements within the associated coarse-scale 
    %element as well as the coarse scale element itself:
    if scale_factor > 1
        temp1 = find(X(:,7) == index_big);
        temp2 = find(X(:,7)==0 & X(:,6)==index_big);
        assoc = [temp1;temp2];
        X_assoc = X(assoc,:);
        X(assoc,:) = [];
        X = [X_assoc;X];
    end
    
    %Truncate to n_max elements
    X = X(1:n_max,:);
    
    %The rest of this function is sim. to above where n_max is not considered
            
    for i = 1:n_max
        %Get next entry in SIGMA12
        X1 = X(i,2);
        Y1 = X(i,3);
        SF1 = X(i,5);
        X2 = indices(4,1);
        Y2 = indices(5,1);
        SF2 = scale_factor;
        
        dX = X2 - X1;
        dY = Y2 - Y1;
        if dY~=0
            dX = sign(dX./dY) * abs(dX);
        end
        dY = abs(dY);
        
        scale_index = (SF1 + SF2) - 1;
        
        %interpolate correlation from previously calculated RHO matrix.
        if SF2 == 1
            rho11 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,5), dX, dY,'*linear');
            rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
            rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
            rho_temp = [rho11,rho12,rho22];
            [r,c] = find(isnan(rho_temp));
            rho_temp(r,c) = 0;
            rho = rho_temp(:,scale_index);
        else
            rho11 = 0;
            rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
            rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
            rho_temp = [rho11,rho12,rho22];
            [r,c] = find(isnan(rho_temp));
            rho_temp(r,c) = 0;
            rho = rho_temp(:,scale_index);
        end 
        
        sigma_temp = [sigma_avg, sigma];
        sigma_n = sigma_temp(:,scale_factor);
        sigma_1 = sigma_temp(:,SF1);
        SIGMA12(1,i) = rho*sigma_n*sigma_1;

        %Get next row in SIGMA22
        X1 = X1 * ones(n_max,1);
        Y1 = Y1 * ones(n_max,1);
        SF1 = SF1 * ones(n_max,1);
        X2 = X(:,2);
        Y2 = X(:,3);
        SF2 = X(:,5);
    
        dX = X2 - X1;
        dY = Y2 - Y1;
        for p = 1:n_max
            if dY(p)~=0
                dX(p) = sign(dX(p)./dY(p)) * abs(dX(p));
            end
            dY(p) = abs(dY(p));
        end

        scale_index = (SF1 + SF2) - 1;                   %vector of SF's

        %interpolate correlation from previously calculated RHO matrix.
        rho11 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,5), dX, dY,'*linear');
        rho12 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,4), dX, dY,'*linear');
        rho22 = interp2(RHO(:,:,1), RHO(:,:,2), RHO(:,:,3), dX, dY,'*linear');
        rho_temp = [rho11,rho12,rho22];
        [r,c] = find(isnan(rho_temp));
        rho_temp(r,c) = 0;
        rho = diag(rho_temp(:,scale_index));
        
        sigma_temp = ones(n_max,1) * [sigma_avg, sigma];
        sigma_1 = ones(n_max,1) * sigma_1;
        sigma_2 = diag(sigma_temp(:,SF2));
        SIGMA22(i,1:n_max) = [rho.*sigma_1.*sigma_2]';
    end
end

%determine conditional parameters
if l <= n_max
    Z_col = USED_MATRIX(:,6);
    modifier = SIGMA12 * (SIGMA22 \ (Z_col-mu));
    mu_prime = mu + modifier;
else
    Z_col = X(:,4);
    modifier = SIGMA12 * (SIGMA22 \ (Z_col(1:n_max)-mu));
    mu_prime = mu + modifier;
end

if scale_factor == 1
    sigma_base = sigma_avg;
else
    sigma_base = sigma;
end

sigma_prime = sqrt(sigma_base^2 - SIGMA12 * (SIGMA22 \ SIGMA12'));

if sigma_prime < 0.0001  %counteracts numerical instabilities at very low conditional variance
    sigma_prime = 0;
end

%end of file