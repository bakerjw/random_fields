%transform.m%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Andy Seifried, Stanford University
%12.27.2009

%This function transforms the simulated fields between the standard normal
%(denoted by 'z') and user-specified (denoted by 'u') distributions, 
%depending on the value of 'trans_flag'.

%Call with: USED_MATRIX and HISTORY (golbal call)
%           FLAG, ds, d_base, a, b, ROTATE, flag, sigma_avg
%Return:    USED_MATRIX and HISTORY (global return)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = transform(trans_flag,ds,d_base,a,b,ROTATE,flag,sigma_avg)

global USED_MATRIX
global HISTORY

%properties defined at fine scale%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Target distribution is lognormal with the following properties:
mu_cu = 100; %kN/m^2
sig_cu = 50;

sig_lncu = sqrt(log(1+(sig_cu/mu_cu)^2));
mu_lncu  = log(mu_cu) - 0.5*sig_lncu^2;

%properties calculated for local averages

var_cu_avg = get_var_avg(ds,d_base,sig_cu,a,b,ROTATE,flag);
sig_cu_avg = sqrt(var_cu_avg);

sig_lncu_avg = sqrt(log(1+(sig_cu_avg/mu_cu)^2));
mu_lncu_avg  = log(mu_cu) - 0.5*sig_lncu_avg^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%perform transformation

if trans_flag == 1 %transform from z to u

    for i = 1:length(USED_MATRIX(:,1))
        if USED_MATRIX(i,7) == 1   % coarse scale
            USED_MATRIX(i,6)=logninv(normcdf(USED_MATRIX(i,6),0,sigma_avg),mu_lncu_avg,sig_lncu_avg);
        else                       % fine scale
            USED_MATRIX(i,6)=logninv(normcdf(USED_MATRIX(i,6),0,1),mu_lncu,sig_lncu);
        end
    end

elseif trans_flag == 2 %transform from u to z
    
    for i = 1:length(USED_MATRIX(:,1))
        if USED_MATRIX(i,7) == 1   % coarse scale
            USED_MATRIX(i,6)=norminv(logncdf(USED_MATRIX(i,6),mu_lncu_avg,sig_lncu_avg),0,sigma_avg);
        else                       % fine scale
            USED_MATRIX(i,6)=norminv(logncdf(USED_MATRIX(i,6),mu_lncu,sig_lncu),0,1);
        end
    end

end

%Also need to transform 'HISTORY':
if trans_flag == 1 %transform from z to u

    for i = 1:length(HISTORY(:,1))
        if HISTORY(i,7) == 1   % coarse scale
            HISTORY(i,6)=logninv(normcdf(HISTORY(i,6),0,sigma_avg),mu_lncu_avg,sig_lncu_avg);
        else                       % fine scale
            HISTORY(i,6)=logninv(normcdf(HISTORY(i,6),0,1),mu_lncu,sig_lncu);
        end
    end

elseif trans_flag == 2 %transform from u to z
    
    for i = 1:length(HISTORY(:,1))
        if HISTORY(i,7) == 1   % coarse scale
            HISTORY(i,6)=norminv(logncdf(HISTORY(i,6),mu_lncu_avg,sig_lncu_avg),0,sigma_avg);
        else                       % fine scale
            HISTORY(i,6)=norminv(logncdf(HISTORY(i,6),mu_lncu,sig_lncu),0,1);
        end
    end
end

%This function transforms the standard normal distribution to the truncated
%exponential distribution used in Andrade, Baker, and Ellison 2007.

%for i = 1:length(USED_MATRIX(:,1))
%    USED_MATRIX(i,6)=(27.5-log(1-normcdf(USED_MATRIX(i,6),0,1)*(1-exp(-5))))/50;
%end

%end of file