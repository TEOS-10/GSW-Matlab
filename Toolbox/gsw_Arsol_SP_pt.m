function Arsol = gsw_Arsol_SP_pt(SP,pt)

% gsw_Arsol_SP_pt                              solubility of Ar in seawater
%==========================================================================
%
% USAGE:  
%  Arsol = gsw_Arsol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the argon, Ar, concentration expected at equilibrium with air 
%  at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including
%  saturated water vapor  This function uses the solubility coefficients
%  as listed in Hamme and Emerson (2004).
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
%  Arsol = solubility of argon                                  [ umol/kg ] 
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

if nargin ~=2
   error('gsw_Arsol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_Arsol_SP_pt: SP and pt must have same dimensions')
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

y = log((298.15 - pt)./(273.15 + pt)); % pt is the temperature in degress C  
                     % on the 1990 International Temperature Scale ITS-90.

% The coefficents below are from Table 4 of Hamme and Emerson (2004)
a0 =  2.79150;
a1 =  3.17609;
a2 =  4.13116;
a3 =  4.90379;
b0 = -6.96233e-3;
b1 = -7.66670e-3;
b2 = -1.16888e-2;

Arsol = exp(a0 + y.*(a1 + y.*(a2 + a3*y)) + x.*(b0 + y.*(b1 + b2*y)));

if transposed
    Arsol = Arsol.';
end

end