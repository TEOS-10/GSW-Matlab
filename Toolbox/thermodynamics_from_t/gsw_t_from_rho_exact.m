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
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
    error('gsw_t_from_rho_exact:  Requires three inputs')
end 

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
end

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

SA(SA<0 | SA>42 | p <-1.5 | p>12000) = NaN;   % SA out of range, set to NaN

rho_40 = gsw_rho_t_exact(SA,40*ones(size(SA)),p);
SA((rho - rho_40) < 0) = NaN;                   % rho too light, set to NaN

if any(isnan(SA(:) + p(:) + rho(:)))
    [I_bad] = find(isnan(SA + p + rho));
    SA(I_bad) = NaN;
end
t_max_rho = gsw_t_maxdensity_exact(SA,p);
max_rho = gsw_rho_t_exact(SA,t_max_rho,p);
v_t_t = gsw_gibbs(0,2,1,SA,t_max_rho,p);
rho_t_t = - v_t_t.*rho.*rho;                          % This is approximate 
deriv_lin = (rho_40 - max_rho)./(40 - t_max_rho);
discrim = deriv_lin.*deriv_lin - 4.*0.5.*rho_t_t.*(max_rho - rho);
discrim(discrim < 0) = 0;
t = t_max_rho + 2.*(max_rho - rho)./(-deriv_lin + sqrt(discrim));

%  Having found this initial value of t, begin the iterative solution
%  for the t_upper part, the solution warmer than the TMD
for Number_of_iterations = 1:3
    v_t = gsw_gibbs(0,1,1,SA,t,p);
    rho_t = -v_t.*rho.*rho;                           % This is approximate
    v_t_t = gsw_gibbs(0,2,1,SA,t,p);
    rho_t_t = -v_t_t.*rho.*rho;                       % This is approximate
       b = rho_t;
       a = 0.5.*rho_t_t;
       c = gsw_rho_t_exact(SA,t,p) - rho;    
       discrim = b.*b - 4.*a.*c;
       discrim(discrim < 0) = 0;
    t_old = t;
    t = t_old + (2.*c)./(-b + sqrt(discrim));
    
end
t_upper = t;
   
%  Now start the t_multiple part, the solution cooler than the TMD
t = 2.*t_max_rho - t_upper;
for Number_of_iterations = 1:6
    v_t = gsw_gibbs(0,1,1,SA,t,p);
    rho_t = -v_t.*rho.*rho;                           % This is approximate
    v_t_t = gsw_gibbs(0,2,1,SA,t,p);
    rho_t_t = -v_t_t.*rho.*rho;                       % This is approximate
       b = rho_t;
       a = 0.5.*rho_t_t;
       c = gsw_rho_t_exact(SA,t,p) - rho;
       discrim = b.*b - 4.*a.*c;
       discrim(discrim < 0) = 0;
    t_old = t;
    t = t_old + (2.*c)./(-b - sqrt(discrim));
                                    % Note the sign change of the sqrt term
end
t_multiple = t;
    
% Set values outside the relevant bounds to NaNs 
t_freezing = gsw_t_freezing(SA,p,0);    % This assumes that the seawater is 
                                        % always unsaturated with air
                                       
t_upper(t_upper < t_freezing) = NaN;
t_upper(t_upper < t_max_rho) = NaN;
t_multiple(t_multiple < t_freezing) = NaN;
t_multiple(t_multiple > t_max_rho) = NaN;

   t = t_upper;

if transposed
    t = t.';
    t_multiple = t_multiple.';
end

end
