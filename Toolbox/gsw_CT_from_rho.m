function [CT,CT_multiple] = gsw_CT_from_rho(rho,SA,p)

% gsw_CT_from_rho                     Conservative Temperature from density
%                                                        (48-term equation)
% =========================================================================
%
% USAGE:
%  [CT,CT_multiple] = gsw_CT_from_rho(rho,SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar), 
%  using the computationally-efficient 48-term expression for density in 
%  terms of SA, CT and p (IOC et al., 2010)
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  CT  =  Conservative Temperature  (ITS-90)                      [ deg C ]
%  CT_multiple  =  Conservative Temperature  (ITS-90)             [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      Conservative Temperatures for a single density.  This programme will
%      output both valid solutions.  To see this second solution the user 
%      must call the programme with two outputs (i.e. [CT,CT_multiple]), if
%      there is only one possible solution and the programme has been 
%      called with two outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
    error('gsw_CT_from_rho:  Requires three inputs')
end %if

[md,nd] = size(rho);
[ms,ns] = size(SA);
[mp,np] = size(p);

if (ms ~= md | ns ~= nd)
    error('gsw_CT_from_rho: rho and SA must have same dimensions')
end

if (mp == 1) & (np == 1)                    % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)               % p is row vector,
    p = p(ones(1,md), :);                   % copy down each column.
elseif (md == mp) & (np == 1)               % p is column vector,
    p = p(:,ones(1,nd));                    % copy across each row.
elseif (nd == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                 % transposed then
    p = p(ones(1,md), :);                   % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_CT_from_rho: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% alpha_limit is the positive value of the thermal expansion coefficient
% which is used at the freezing temperature to distinguish between
% I_salty and I_fresh.
alpha_limit = 1e-5;

% rec_half_rho_TT is a constant representing the reciprocal of half the
% second derivative of density with respect to temperature near the
% temperature of maximum density.
rec_half_rho_TT = -110.0;

CT = nan(size(SA));
CT_multiple = CT;

% SA out of range, set to NaN.
SA(SA<0 | SA>42 | p <-1.5 | p>12000) = NaN;

rho_40 = gsw_rho_CT(SA,40*ones(size(SA)),p);
% rho too light, set to NaN.
SA((rho - rho_40) < 0) = NaN;

CT_max_rho = gsw_CT_maxdensity(SA,p);
rho_max = gsw_rho(SA,CT_max_rho,p);
rho_extreme = rho_max;
CT_freezing = gsw_CT_freezing(SA,p); % this assumes that the seawater is always saturated with air
rho_freezing = gsw_rho(SA,CT_freezing,p);
% reset the extreme values
rho_extreme((CT_freezing - CT_max_rho) > 0) = rho_freezing((CT_freezing - CT_max_rho) > 0);

% set SA values to NaN for the rho's that are too dense.
SA(rho > rho_extreme) = NaN;

if any(isnan(SA + p + rho))
    [I_bad] = find(isnan(SA + p + rho));
    SA(I_bad) = NaN;
end

alpha_freezing = gsw_alpha(SA,CT_freezing,p);

if any(alpha_freezing > alpha_limit)
    [I_salty] = find(alpha_freezing > alpha_limit);
    CT_diff = 40*ones(size(I_salty)) - CT_freezing(I_salty);
    
    top = rho_40(I_salty) - rho_freezing(I_salty) ...
           + rho_freezing(I_salty).*alpha_freezing(I_salty).*CT_diff;
    a = top./(CT_diff.*CT_diff);
    b = - rho_freezing(I_salty).*alpha_freezing(I_salty);
    c = rho_freezing(I_salty) - rho(I_salty);
    sqrt_disc = sqrt(b.*b - 4*a.*c);
    % the value of t(I_salty) here is the initial guess at CT in the range 
    % of I_salty.
    CT(I_salty) = CT_freezing(I_salty) + 0.5*(-b - sqrt_disc)./a;
end

if any(alpha_freezing <= alpha_limit)
    [I_fresh] = find(alpha_freezing <= alpha_limit);

    CT_diff = 40*ones(size(I_fresh)) - CT_max_rho(I_fresh);
    factor = (rho_max(I_fresh) - rho(I_fresh))./ ...
               (rho_max(I_fresh) - rho_40(I_fresh));
    delta_CT = CT_diff.*sqrt(factor);
    
    if any(delta_CT > 5)
        [I_fresh_NR] = find(delta_CT > 5);
        CT(I_fresh(I_fresh_NR)) = CT_max_rho(I_fresh(I_fresh_NR)) + delta_CT(I_fresh_NR);
    end
     
    if any(delta_CT <= 5)
        [I_quad] = find(delta_CT <= 5);
        CT_a = nan(size(SA));
        % set the initial value of the quadratic solution routes.
        CT_a(I_fresh(I_quad)) = CT_max_rho(I_fresh(I_quad)) + ...
            sqrt(rec_half_rho_TT*(rho(I_fresh(I_quad)) - rho_max(I_fresh(I_quad))));       
        for Number_of_iterations = 1:7
            CT_old = CT_a;
            rho_old = gsw_rho(SA,CT_old,p);
            factorqa = (rho_max - rho)./(rho_max - rho_old);
            CT_a = CT_max_rho + (CT_old - CT_max_rho).*sqrt(factorqa);
        end
        
        CT_a(CT_freezing - CT_a < 0) = NaN;
        
        CT_b = nan(size(SA));
        % set the initial value of the quadratic solution roots.
        CT_b(I_fresh(I_quad)) = CT_max_rho(I_fresh(I_quad)) - ...
            sqrt(rec_half_rho_TT*(rho(I_fresh(I_quad)) - rho_max(I_fresh(I_quad))));    
        for Number_of_iterations = 1:7
            CT_old = CT_b;
            rho_old = gsw_rho(SA,CT_old,p);
            factorqb = (rho_max - rho)./(rho_max - rho_old);
            CT_b = CT_max_rho + (CT_old - CT_max_rho).*sqrt(factorqb);
        end
% After seven iterations of this quadratic iterative procedure,
% the error in rho is no larger than 4.6x10^-13 kg/m^3.
        CT_b(CT_freezing - CT_b < 0) = NaN;
    end
end

% begin the modified Newton-Raphson iterative method, which will only
% operate on non-NaN CT data.

v_lab = ones(size(rho))./rho;
v_CT = gsw_specvol(SA,CT,p).*gsw_alpha(SA,CT,p);

for Number_of_iterations = 1:3
    CT_old = CT;
    delta_v = gsw_specvol(SA,CT_old,p) - v_lab;
    CT = CT_old - delta_v./v_CT ; % this is half way through the modified N-R method
    CT_mean = 0.5*(CT + CT_old);
    v_CT = gsw_specvol(SA,CT_mean,p).*gsw_alpha(SA,CT_mean,p);
    CT = CT_old - delta_v./v_CT ;
end

if exist('t_a','var')
    CT(~isnan(CT_a)) = CT_a(~isnan(CT_a));
end
if exist('t_b','var')
    CT_multiple(~isnan(CT_b)) = CT_b(~isnan(CT_b));
end
% After three iterations of this modified Newton-Raphson iteration,
% the error in rho is no larger than 1.6x10^-12 kg/m^3.

if transposed
    CT = CT.';
    CT_multiple = CT_multiple.';
end

end
