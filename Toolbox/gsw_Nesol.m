function Nesol = gsw_Nesol(SA,CT,p,long,lat)

% gsw_Nesol                                    solubility of Ne in seawater
%==========================================================================
%
% USAGE:  
%  Nesol = gsw_Nesol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the Neon, Ne, concentration expected at equilibrium with air 
%  at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor.  This function uses the solubility coefficients
%  as listed in Hamme and Emerson (2004).
%
%  Note that this algorithm has not been approved by IOC and is not work 
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:  
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions. p, lat and long may have 
%  dimensions 1x1 or Mx1 or 1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  Nesol = solubility of neon in nano-moles per kg              [ nmol/kg ] 
% 
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Hamme, R., and S. Emerson, 2004: The solubility of neon, nitrogen and
%   argon in distilled water and seawater. Deep-Sea Research, 51, 
%   1517-1528.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if nargin ~= 5
   error('gsw_Nesol: Requires five inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_Nesol: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)                 % p scalar - fill to size of SA
    p = p*ones(ms,ns);
elseif (ns == np) & (mp == 1)            % p is row vector,
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (np == 1)            % p is column vector,
    p = p(:,ones(1,ns));                 % copy across each row.
elseif (ns == mp) & (np == 1)            % p is a transposed row vector,
    p = p.';                              % transpose, then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Nesol: Inputs array dimensions arguments do not agree')
end %if

[mla,nla] = size(lat);

if (mla == 1) & (nla == 1)             % lat is a scalar - fill to size of SA
    lat = lat*ones(ms,ns);
elseif (ns == nla) & (mla == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (nla == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));           % copy across each row.
elseif (ns == mla) & (nla == 1)        % lat is a transposed row vector,
    lat = lat.';                        % transpose, then
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mla) & (ns == nla)
    % ok
else
    error('gsw_Nesol: Inputs array dimensions arguments do not agree')
end %if

[mlo,nlo] = size(long);
long(long < 0) = long(long < 0) + 360; 

if (mlo == 1) & (nlo == 1)            % long is a scalar - fill to size of SA
    long = long*ones(ms,ns);
elseif (ns == nlo) & (mlo == 1)       % long is a row vector,
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (nlo == 1)       % long is a column vector,
    long = long(:,ones(1,ns));        % copy across each row.
elseif (ns == mlo) & (nlo == 1)       % long is a transposed row vector,
    long = long.';                     % transpose, then
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mlo) & (ns == nlo)
    % ok
else
    error('gsw_Nesol: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    lat = lat.';
    long = long.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SP = gsw_SP_from_SA(SA,p,long,lat);

x = SP;        % Note that salinity argument is Practical Salinity, this is
             % beacuse the major ionic components of seawater related to Cl  
          % are what affect the solubility of non-electrolytes in seawater.   
                 
pt = gsw_pt_from_CT(SA,CT); % pt is potential temperature referenced to
                            % the sea surface.
y = log((298.15 - pt)./(gsw_T0 + pt)); % pt is the temperature in degress C  
                     % on the 1990 International Temperature Scale ITS-90.

% The coefficents below are from Table 4 of Hamme and Emerson (2004)
a0 =  2.18156;
a1 =  1.29108;
a2 =  2.12504;
b0 = -5.94737e-3;
b1 = -5.13896e-3;

Nesol = exp(a0 + y.*(a1 + a2*y) + x.*(b0 + b1.*y));

if transposed
    Nesol = Nesol.';
end

end