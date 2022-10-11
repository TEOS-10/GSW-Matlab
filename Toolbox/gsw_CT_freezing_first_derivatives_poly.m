function [CTfreezing_SA, CTfreezing_P] = gsw_CT_freezing_first_derivatives_poly(SA,p,saturation_fraction)

% gsw_CT_freezing_first_derivatives_poly               first derivatives of
%                 Conservative Temperature at which seawater freezes (poly)
%==========================================================================
%
% USAGE:
%  [CTfreezing_SA, CTfreezing_P] = gsw_CT_freezing_first_derivatives_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the Conservative Temperature at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa) of the comptationally efficient polynomial fit of the
%  freezing temperature (McDougall et al., 2014).
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
%  CTfreezing_SA = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to Absolute Salinity at
%                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  CTfreezing_P  = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to pressure (in Pa) at
%                  fixed Absolute Salinity                         [ K/Pa ]
%
% AUTHOR: 
%  Trevor McDougall, Paul Barker  [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
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
   error('gsw_CT_freezing_first_derivatives_poly:  Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_CT_freezing_first_derivatives_poly: saturation_fraction MUST be between zero and one.')
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
    error('gsw_CT_freezing_first_derivatives_poly: Inputs array dimensions arguments do not agree')
end %if

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
    error('gsw_CT_freezing_first_derivatives_poly: Inputs array dimensions arguments do not agree')
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

SA(SA < 0) = 0; % This line ensures that SA is not negative

% c0  =  0.017947064327968736;
c1 =  -6.076099099929818;
c2 =   4.883198653547851;
c3 =  -11.88081601230542;
c4 =   13.34658511480257;
c5 =  -8.722761043208607;
c6 =   2.082038908808201;
%
c7 =  -7.389420998107497;
c8 =  -2.110913185058476;
c9 =   0.2295491578006229; 
%  
c10 = -0.9891538123307282;
c11 = -0.08987150128406496;
c12 =  0.3831132432071728;
c13 =  1.054318231187074;
c14 =  1.065556599652796;
c15 = -0.7997496801694032;
c16 =  0.3850133554097069;
c17 = -2.078616693017569;
c18 =  0.8756340772729538;
c19 = -2.079022768390933;
c20 =  1.596435439942262;
c21 =  0.1338002171109174;
c22 =  1.242891021876471;

SA_r = SA.*1e-2;
x = sqrt(SA_r);
p_r = p.*1e-4;

% a = 0.014289763856964;             % Note that a = 0.502500117621/35.16504.
% b = 0.057000649899720;
% Note that -0.018994561378548 = -a -a*b - 2.4*b/gsw_SSO
% and, 4.632588654871302e-05 = 2*a*b./gsw_SSO

CTfreezing_SA = (c1 + x.*(1.5*c2  + x.*(2*c3  + x.*(2.5*c4  + x.*(3*c5  + 3.5*c6*x)))) ...
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22*x)))) ... 
       + p_r.*(c12 + x.*(1.5*c14 + x.*(2*c17 + 2.5*c20*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21*x))))).*1e-2...
        - saturation_fraction.*(1e-3).*(-0.018994561378548 - SA.*4.632588654871302e-05);

CTfreezing_P = (c7 + SA_r.*(c10   + x.*(c11   + x.*(c13 + x.*(c16 + x.*(c19 + c22*x))))) ...
     + p_r.*(2*c8 + SA_r.*(2*c12 + x.*(2*c14 + x.*(2*c17 + 2*c20*x))) ...
     + p_r.*(3*c9 + SA_r.*(3*c15 + x.*(3*c18 + 3*c21*x))))).*1e-8;

% set any values that are out of range to be NaN. 
CTfreezing_SA(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

CTfreezing_P(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    CTfreezing_SA = CTfreezing_SA.';
    CTfreezing_P = CTfreezing_P.';
end

end