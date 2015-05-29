function Hesol = gsw_Hesol_SP_pt(SP,pt)

% gsw_Hesol_SP_pt                              solubility of He in seawater
%==========================================================================
%
% USAGE:  
%  Hesol = gsw_Hesol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the helium concentration expected at equilibrium with air at 
%  an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) including 
%  saturated water vapor.  This function uses the solubility coefficients
%  as listed in Weiss (1971).
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
%  Hesol = solubility of helium in micro-moles per kg           [ umol/kg ] 
% 
% AUTHOR:  Roberta Hamme, Paul Barker and Trevor McDougall
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Dymond and Smith, 1980: The virial coefficients of pure gases and
%   mixtures. Clarendon Press, Oxford.
%
%  Weiss, R.F., 1971: Solubility of Helium and Neon in Water and Seawater.
%   J. Chem. and Engineer. Data, 16, 235-241.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if nargin ~=2
   error('gsw_Hesol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_Hesol_SP_pt: SP and pt must have same dimensions')
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

pt68 = pt.*1.00024; % pt68 is the potential temperature in degress C on 
              % the 1968 International Practical Temperature Scale IPTS-68.
y = pt68 + gsw_T0;
y_100 = y.*1e-2;

% The coefficents below are from Table 3 of Weiss (1971)
a1 = -167.2178;
a2 =  216.3442;
a3 =  139.2032;
a4 = -22.6202;
b1 = -0.044781;
b2 =  0.023541;
b3 = -0.0034266;

Hesol = exp(a1 + a2*100./y + a3*log(y_100) + a4*y_100 ...
           + x.*(b1 + y_100.*(b2 + b3*y_100)));

He_ml2umol = 4.455817671505537e1; % mL/kg to umol/kg for He (1/22.44257e-3) 
                             %Molar volume at STP (Dymond and Smith, 1980).
Hesol = Hesol.*He_ml2umol;

if transposed
    Hesol = Hesol.';
end

end