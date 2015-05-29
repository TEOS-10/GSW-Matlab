function [t,t_multiple] = gsw_t_from_rho_exact(rho,SA,p)

% gsw_t_from_rho_exact                     in-situ temperature from density
% =========================================================================
%
% USAGE:
%  [t,t_multiple] = gsw_t_from_rho_exact(rho,SA,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg m^-3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  t  =  in-situ temperature  (ITS-90)                            [ deg C ]
%  t_multiple  =  in-situ temperature  (ITS-90)                   [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both 
%      valid solutions.  To see this second solution the user must call the 
%      programme with two outputs (i.e. [t,t_multiple]), if there is only 
%      one possible solution and the programme has been called with two 
%      outputs the second variable will be set to NaN.
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
%  McDougall, T.J. and S.J. Wotherspoon, 2012: A simple modification of 
%   Newton’s method to achieve convergence of order "1 + sqrt(2)".
%   Submitted to Applied Mathematics and Computation.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
    error('gsw_t_from_rho_exact:  Requires three inputs')
end %if

[md,nd] = size(rho);
[ms,ns] = size(SA);
[mp,np] = size(p);

if (ms ~= md | ns ~= nd)
    error('gsw_t_from_rho_exact: rho and SA must have same dimensions')
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
    error('gsw_t_from_rho_exact: Inputs array dimensions arguments do not agree')
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

t = nan(size(SA));
t_multiple = t;

SA(SA<0 | SA>42 | p <-1.5 | p>12000) = NaN;

rho_40 = gsw_rho_t_exact(SA,40*ones(size(SA)),p);
    
SA((rho - rho_40) < 0) = NaN;

t_max_rho = gsw_t_maxdensity_exact(SA,p);
rho_max = gsw_rho_t_exact(SA,t_max_rho,p);
rho_extreme = rho_max;
t_freezing = gsw_t_freezing(SA,p); % this assumes that the seawater is always saturated with air
rho_freezing = gsw_rho_t_exact(SA,t_freezing,p);

%set rhos greater than those at the freexing point to be equal to the freezing point.
rho_extreme((t_freezing - t_max_rho) > 0) = rho_freezing((t_freezing - t_max_rho) > 0);
% set rho that are greater than the extreme limits to NaN.
SA(rho > rho_extreme) = NaN;

SA(isnan(SA + p + rho)) = NaN;

alpha_freezing = gsw_alpha_wrt_t_exact(SA,t_freezing,p);

if any(alpha_freezing > alpha_limit)
    [I_salty] = find(alpha_freezing > alpha_limit);

    t_diff = 40*ones(size(I_salty)) - t_freezing(I_salty);
    top = rho_40(I_salty) - rho_freezing(I_salty) ...
           + rho_freezing(I_salty).*alpha_freezing(I_salty).*t_diff;
    a = top./(t_diff.*t_diff);
    b = - rho_freezing(I_salty).*alpha_freezing(I_salty);
    c = rho_freezing(I_salty) - rho(I_salty);
    sqrt_disc = sqrt(b.*b - 4*a.*c);
    % the value of t(I_salty) here is the initial guess at t in the range of
    % I_salty.
    t(I_salty) = t_freezing(I_salty) + 0.5*(-b - sqrt_disc)./a;
end

if any(alpha_freezing <= alpha_limit)
    [I_fresh] = find(alpha_freezing <= alpha_limit);

    t_diff = 40*ones(size(I_fresh)) - t_max_rho(I_fresh);
    factor = (rho_max(I_fresh) - rho(I_fresh))./ ...
               (rho_max(I_fresh) - rho_40(I_fresh));
    delta_t = t_diff.*sqrt(factor);
    
    if any(delta_t > 5)
        [I_fresh_NR] = find(delta_t > 5);
        t(I_fresh(I_fresh_NR)) = t_max_rho(I_fresh(I_fresh_NR)) + delta_t(I_fresh_NR);
    end
    
    if any(delta_t <= 5)
       [I_quad] = find(delta_t <= 5);
        
        t_a = nan(size(SA));
        % set the initial value of the quadratic solution roots.
        t_a(I_fresh(I_quad)) = t_max_rho(I_fresh(I_quad)) + ...
            sqrt(rec_half_rho_TT*(rho(I_fresh(I_quad)) - rho_max(I_fresh(I_quad))));       
        for Number_of_iterations = 1:7
            t_old = t_a;
            rho_old = gsw_rho_t_exact(SA,t_old,p);
            factorqa = (rho_max - rho)./(rho_max - rho_old);
            t_a = t_max_rho + (t_old - t_max_rho).*sqrt(factorqa);
        end        
% set temperatures that are less than the freezing point to NaN 
        t_a(t_freezing - t_a < 0) = NaN;
        
        t_b = nan(size(SA));
        % set the initial value of the quadratic solution routes.
        t_b(I_fresh(I_quad)) = t_max_rho(I_fresh(I_quad)) - ...
            sqrt(rec_half_rho_TT*(rho(I_fresh(I_quad)) - rho_max(I_fresh(I_quad))));    
        for Number_of_iterations = 1:7
            t_old = t_b;
            rho_old = gsw_rho_t_exact(SA,t_old,p);
            factorqb = (rho_max - rho)./(rho_max - rho_old);
            t_b = t_max_rho + (t_old - t_max_rho).*sqrt(factorqb);
        end
% After seven iterations of this quadratic iterative procedure,
% the error in rho is no larger than 4.6x10^-13 kg/m^3.
% set temperatures that are less than the freezing point to NaN 
        t_b(t_freezing - t_b < 0) = NaN;
    end
end

% begin the modified Newton-Raphson iterative method (McDougall and 
% Wotherspoon, 2012), which will only operate on non-NaN t data.

v_lab = ones(size(rho))./rho;
v_t = gsw_gibbs(0,1,1,SA,t,p);
for Number_of_iterations = 1:4
    t_old = t;
    delta_v = gsw_gibbs(0,0,1,SA,t_old,p) - v_lab;
    t = t_old - delta_v./v_t ; % this is half way through the modified N-R method
    t_mean = 0.5*(t + t_old);
    v_t = gsw_gibbs(0,1,1,SA,t_mean,p);
    t = t_old - delta_v./v_t ;
end

if exist('t_a','var')
    t(~isnan(t_a)) = t_a(~isnan(t_a));
end
if exist('t_b','var')
    t_multiple(~isnan(t_b)) = t_b(~isnan(t_b));
end
% After three iterations of this modified Newton-Raphson iteration,
% the error in rho is no larger than 4.6x10^-13 kg/m^3.

if transposed
    t = t.';
    t_multiple = t_multiple.';
end

end
