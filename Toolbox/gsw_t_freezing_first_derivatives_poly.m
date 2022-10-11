function [tfreezing_SA, tfreezing_P] = gsw_t_freezing_first_derivatives_poly(SA,p,saturation_fraction)

% gsw_t_freezing_first_derivatives_poly            first derivatives of the  
%                             in-situ temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  [tfreezing_SA, tfreezing_P] = ...
%           gsw_t_freezing_first_derivatives_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the frist derivatives of the in-situ temperature at which 
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
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_t_freezing_first_derivatives_poly: Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_t_freezing_first_derivatives_poly: saturation fraction MUST be between zero and one.')
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
    error('gsw_t_freezing_first_derivatives_poly: Inputs array dimensions arguments do not agree')
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
    error('gsw_t_freezing_first_derivatives_poly: Inputs array dimensions arguments do not agree')
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

%c0 = 0.002519;

c1 = -5.946302841607319;
c2 =  4.136051661346983;
c3 = -1.115150523403847e1;
c4 =  1.476878746184548e1;
c5 = -1.088873263630961e1;
c6 =  2.961018839640730;
    
c7 = -7.433320943962606;
c8 = -1.561578562479883;
c9 =  4.073774363480365e-2;

c10 =  1.158414435887717e-2;
c11 = -4.122639292422863e-1;
c12 = -1.123186915628260e-1;
c13 =  5.715012685553502e-1;
c14 =  2.021682115652684e-1;
c15 =  4.140574258089767e-2;
c16 = -6.034228641903586e-1;
c17 = -1.205825928146808e-2;
c18 = -2.812172968619369e-1;
c19 =  1.877244474023750e-2;
c20 = -1.204395563789007e-1;
c21 =  2.349147739749606e-1;
c22 =  2.748444541144219e-3;

SA_r = SA.*1e-2;
x = sqrt(SA_r);
p_r = p.*1e-4;

tfreezing_SA = (c1 + x.*(1.5*c2 + x.*(2*c3 + x.*(2.5*c4 + x.*(3*c5 + 3.5*c6.*x))))  ... 
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22.*x)))) ...
       + p_r.*(c12 + x.*(1.5*c14 + x.*(2*c17 + 2.5*c20.*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21.*x))))).*1e-2 ...
       + 1.421866717626370e-05.*saturation_fraction;
   
tfreezing_P = (c7 + SA_r.*(c10 + x.*(c11 + x.*(c13 + x.*(c16 + x.*(c19 + c22*x))))) ...
    + p_r.*(2*c8 + SA_r.*(2*c12 + x.*(2*c14 + x.*(2*c17+ 2*c20.*x))) ...
    + p_r.*(3*c9 + SA_r.*(3*c15 + x.*(3*c18 + 3*c21.*x)))))*1e-8;

if transposed
    tfreezing_SA = tfreezing_SA.';
    tfreezing_P = tfreezing_P.';
end

end