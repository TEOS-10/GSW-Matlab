function Krsol = gsw_Krsol(SA,CT,p,long,lat)

% gsw_Krsol                                    solubility of Kr in seawater
%==========================================================================
%
% USAGE:  
%  Krsol = gsw_Krsol(SA,CT,p,long,lat)
%
% DESCRIPTION:
%  Calculates the krypton, Kr, concentration expected at equilibrium with  
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) 
%  including saturated water vapor.  This function uses the solubility 
%  coefficients derived from the data of Weiss and Kyser (1978).
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
%  Krsol = solubility of krypton in micro-moles per kg          [ umol/kg ] 
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
%  Weiss, R.F., and T.K. Kyser, 1978: Solubility of Krypton in Water and 
%   Seawater. J. Chem. Thermodynamics, 23, 69-72.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if nargin ~=5
   error('gsw_Krsol: Requires five inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_Krsol: SA and CT must have same dimensions')
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
    error('gsw_Krsol: Inputs array dimensions arguments do not agree')
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
    error('gsw_Krsol: Inputs array dimensions arguments do not agree')
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
    error('gsw_Krsol: Inputs array dimensions arguments do not agree')
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
pt68 = pt.*1.00024; % pt68 is the potential temperature in degress C on 
              % the 1968 International Practical Temperature Scale IPTS-68.
y = pt68 + gsw_T0;
y_100 = y.*1e-2;

% Table 2 (Weiss and Kyser, 1978)
a1 = -112.6840;
a2 =  153.5817;
a3 =  74.4690;
a4 = -10.0189;
b1 = -0.011213;
b2 = -0.001844;
b3 =  0.0011201;

Krsol = exp(a1 + a2*100./y + a3*log(y_100) + a4*y_100 ...
          + x.*(b1 + y_100.*(b2 + b3*y_100)));

Kr_ml2umol = 4.474052731185490e1;  % mL/kg to umol/kg for Kr (1/22.3511e-3)
                             %Molar volume at STP (Dymond and Smith, 1980).
Krsol = Krsol.*Kr_ml2umol;

if transposed
    Krsol = Krsol.';
end

end