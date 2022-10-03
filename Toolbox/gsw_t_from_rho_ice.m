function t_ice = gsw_t_from_rho_ice(rho_ice,p)

% gsw_t_from_rho_ice                        temperature from density of ice
% =========================================================================
%
% USAGE:
%   t_ice = gsw_t_from_rho_ice(rho_ice,p)
% 
% DESCRIPTION:
%  Calculates the in-situ temperature of ice, for given values of its 
%  density and sea pressure (in dbar). 
%
% INPUT:
%  rho_ice =  density of ice                                     [ kg/m^3 ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho_ice is MxN.
%
% OUTPUT:
%  t_ice  =  in-situ temperature of ice                           [ deg C ]
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
%    See section 2.5 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==2)
   error('gsw_t_from_rho_ice:  Requires two inputs')
end %if

[md,nd] = size(rho_ice);
[mp,np] = size(p);

if (mp == 1) & (np == 1)               % p scalar - fill to size of rho
    p = p*ones(md,nd);
elseif (nd == np) & (mp == 1)          % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)          % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                            % transposed then
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_t_from_rho_ice: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho_ice = rho_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

specvol_ice = 1./rho_ice;

v_173_15 = gsw_specvol_ice(-100*ones(size(specvol_ice)),p);
t_freezing = gsw_t_freezing(zeros(size(specvol_ice)),p,0);
v_freezing = gsw_specvol_ice(t_freezing,p);

t_ice = nan(size(rho_ice));
v_t_ice = t_ice;
if any(specvol_ice(:) < v_173_15(:)) 
    v_93_15= gsw_specvol_ice(-180*ones(size(specvol_ice)),p);
    
    [Icolder] = find(specvol_ice < v_93_15);
    v_63_15= gsw_specvol_ice(-210*ones(size(specvol_ice)),p);
    t_ice(Icolder) = -180 + 30*(specvol_ice(Icolder) - v_93_15(Icolder))./(v_93_15(Icolder) - v_63_15(Icolder));  % initial estimate.
    v_t_ice(Icolder) = (v_63_15(Icolder) - v_93_15(Icolder))./(-30); %initial estimate of v_t_ice, the t derivative of v

    [Icold] = find(specvol_ice >= v_93_15 & specvol_ice < v_173_15);
    t_ice(Icold) = -100 + 80*(specvol_ice(Icold) - v_173_15(Icold))./(v_173_15(Icold) - v_93_15(Icold));  % initial estimate.
    v_t_ice(Icold) = (v_93_15(Icold) - v_173_15(Icold))./(-80); %initial estimate of v_t_ice, the t derivative of v
    
    [Istd] = find(specvol_ice >= v_173_15);
    t_ice(Istd) =  (100 - t_freezing(Istd)).*(specvol_ice(Istd) - v_freezing(Istd))./(v_freezing(Istd) - v_173_15(Istd));  % initial estimate.
    v_t_ice(Istd) = (v_freezing(Istd) - v_173_15(Istd))./(-100 - t_freezing(Istd)); %initial estimate of v_t_ice, the t derivative of v
else
    t_ice = (100 - t_freezing).*(specvol_ice - v_freezing)./(v_freezing - v_173_15);  % the initial estimate of t_ice
    v_t_ice = (v_freezing - v_173_15)./(-100 - t_freezing); %initial estimate of v_t_ice, the t derivative of v
end

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

for Number_of_iterations = 1:3 
    t_ice_old = t_ice;
    delta_v_ice = gsw_specvol_ice(t_ice_old,p) - specvol_ice;
    t_ice = t_ice_old - delta_v_ice./v_t_ice ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)
    t_ice_mean = 0.5*(t_ice + t_ice_old);
    v_t_ice = gsw_gibbs_ice(1,1,t_ice_mean,p);
    t_ice = t_ice_old - delta_v_ice./v_t_ice; 
end

% After 3 iterations of this modified Newton-Raphson iteration,
% the error in t_ice is no larger than 2.3x10^-12 deg C, which 
% is machine precision for this calculation. 
 
if transposed
    t_ice = t_ice.';
end

end
