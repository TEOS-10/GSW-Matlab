function N2Osol = gsw_N2Osol_SP_pt(SP,pt)

% gsw_N2Osol_SP_pt                            solubility of N2O in seawater
%==========================================================================
%
% USAGE:  
%  N2Osol = gsw_N2Osol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the nitrous oxide, N2O, concentration expected at equilibrium  
%  with air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) 
%  including saturated water vapor  This function uses the solubility 
%  coefficients as listed in Weiss and Price (1980).
%
%  Note that this algorithm has not been approved by IOC and is not work 
%  from SCOR/IAPSO Working Group 127. It is included in the GSW
%  Oceanographic Toolbox as it seems to be oceanographic best practice.
%
% INPUT:  
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
%
% OUTPUT:
%  N2Osol = solubility of nitrous oxide                           [ mol/L ] 
% 
% AUTHOR:  Rich Pawlowicz, Paul Barker and Trevor McDougall
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
%  Weiss, R.F., and B.A. Price, 1980: Nitrous oxide solubility in water and
%   seawater. Mar. Chem., 8, 347-359.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if nargin ~=2
   error('gsw_N2Osol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_N2Osol_SP_pt: SP and pt must have same dimensions')
end

if ms == 1
    SP = SP';
    pt = pt';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

x = SP;        % Note that salinity argument is Practical Salinity, this is
             % beacuse the major ionic components of seawater related to Cl  
          % are what affect the solubility of non-electrolytes in seawater.   

pt68 = pt.*1.00024; % pt68 is the potential temperature in degress C on 
              % the 1968 International Practical Temperature Scale IPTS-68.
y = pt68 + gsw_T0;
y_100 = y.*1e-2;

% The coefficents below are from Table 2 of Weiss and Price (1980)
a0 = -165.8806;
a1 =  222.8743;
a2 =  92.0792;
a3 = -1.48425;
b1 = -0.056235;
b2 =  0.031619;
b3 = -0.0048472;

 
m0 = 24.4543;
m1 = 67.4509;
m2 = 4.8489;
m3 = 0.000544;

ph2odP = exp(m0 - m1*100./y - m2*log(y_100) - m3*x); % Moist air correction at 1 atm.

N2Osol = (exp(a0 + a1*100./y + a2*log(y_100) + a3*y_100.^2 ...
           + x.*(b1 + y_100.*(b2 + b3*y_100))))./(1-ph2odP);

if transposed
    N2Osol = N2Osol.';
end

end