function [tfreezing_SA, tfreezing_P] = gsw_t_freezing_first_derivatives(SA,p,saturation_fraction)

% gsw_t_freezing_first_derivatives         first derivatives of the in-situ 
%                                     temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  [tfreezing_SA, tfreezing_P] = ...
%                gsw_t_freezing_first_derivatives(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the in-situ temperature at which 
%  seawater freezes with respect to Absolute Salinity SA and pressure P (in
%  Pa).  These expressions come from differentiating the expression that
%  defines the freezing temperature, namely the equality between the 
%  chemical potentials of water in seawater and in ice.  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  tfreezing_SA = the derivative of the in-situ freezing temperature 
%                 (ITS-90) with respect to Absolute Salinity at fixed    
%                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]               
%
%  tfreezing_P  = the derivative of the in-situ freezing temperature  
%                 (ITS-90) with respect to pressure (in Pa) at fixed  
%                 Absolute Salinity                                [ K/Pa ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_t_freezing_first_derivatives: Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_t_freezing_first_derivatives: saturation fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_t_freezing_first_derivatives: Inputs array dimensions arguments do not agree')
end 

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_t_freezing_first_derivatives: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

g_per_kg = 1000; % This is a conversion factor, being the number of g per kg.  

tf = gsw_t_freezing(SA,p,saturation_fraction) ;
rec_denom = 1./(g_per_kg.*gsw_t_deriv_chem_potential_water_t_exact(SA,tf,p) ...
                  + gsw_entropy_ice(tf,p));

tfreezing_SA = gsw_dilution_coefficient_t_exact(SA,tf,p).*rec_denom ...
               + saturation_fraction.*(1e-3)./70.33008;
tfreezing_P = -(gsw_specvol_t_exact(SA,tf,p) - SA.*gsw_gibbs(1,0,1,SA,tf,p) ...
               - gsw_specvol_ice(tf,p)).*rec_denom;

if transposed
    tfreezing_SA = tfreezing_SA.';
    tfreezing_P = tfreezing_P.';
end

end