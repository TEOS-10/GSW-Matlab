function N2sol = gsw_N2sol_SP_pt(SP,pt)

% gsw_N2sol_SP_pt                              solubility of N2 in seawater
%==========================================================================
%
% USAGE:  
%  N2sol = gsw_N2sol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the nitrogen, N2, concentration expected at equilibrium with  
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) 
%  including saturated water vapor.  This function uses the solubility 
%  coefficients as listed in Hamme and Emerson (2004).
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
%  N2sol = solubility of nitrogen in micro-moles per kg         [ umol/kg ] 
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
   error('gsw_N2sol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_N2sol_SP_pt: SP and pt must have same dimensions')
end

if ms == 1
    SP = SP.';
    pt = pt.';
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

y = log((298.15 - pt)./(gsw_T0 + pt)); % pt is the temperature in degress C  
                     % on the 1990 International Temperature Scale ITS-90.

% The coefficents below are from Table 4 of Hamme and Emerson (2004)
a0 = 6.42931;
a1 = 2.92704;
a2 = 4.32531;
a3 = 4.69149;
b0 = -7.44129e-3;
b1 = -8.02566e-3;
b2 = -1.46775e-2;

N2sol = exp(a0 + y.*(a1 + y.*(a2 + a3*y)) + x.*(b0 + y.*(b1 + b2*y)));

if transposed
    N2sol = N2sol.';
end

end