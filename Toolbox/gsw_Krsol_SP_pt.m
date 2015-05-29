function Krsol = gsw_Krsol_SP_pt(SP,pt)

% gsw_Krsol_SP_pt                              solubility of Kr in seawater
%==========================================================================
%
% USAGE:  
%  Krsol = gsw_Krsol_SP_pt(SP,pt)
%
% DESCRIPTION:
%  Calculates the krypton, Kr, concentration expected at equilibrium with  
%  air at an Absolute Pressure of 101325 Pa (sea pressure of 0 dbar) 
%  including saturated water vapor.  This function uses the solubility 
%  coefficients derived from the data of Weiss (1971).
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
%  Krsol = solubility of krypton in micro-moles per kg          [ umol/kg ] 
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
%  Weiss, R.F. and T.K. Kyser, 1978: Solubility of Krypton in Water and 
%   Seawater. J. Chem. Thermodynamics, 23, 69-72.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if nargin ~=2
   error('gsw_Krsol_SP_pt: Requires two inputs')
end %if

[ms,ns] = size(SP);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_Krsol_SP_pt: SP and pt must have same dimensions')
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
y = pt68 +gsw_T0;
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