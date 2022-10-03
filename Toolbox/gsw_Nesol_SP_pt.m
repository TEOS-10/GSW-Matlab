function Nesol = gsw_Nesol_SP_pt(SP,pt)

% gsw_Nesol_SP_pt                              solubility of Ne in seawater
%==========================================================================
%
% USAGE:  
%  Nesol = gsw_Nesol_SP_pt(SP,pt)
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
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  pt  =  potential temperature (ITS-90) referenced               [ deg C ]
%         to one standard atmosphere (0 dbar).
%
%  SP & pt need to have the same dimensions.
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

if nargin ~=2
   error('gsw_Nesol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_Nesol_SP_pt: SP and pt must have same dimensions')
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