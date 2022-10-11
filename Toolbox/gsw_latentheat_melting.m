function latentheat_melting = gsw_latentheat_melting(SA,p)

% gsw_latentheat_melting                             latent heat of melting 
%==========================================================================
%
% USAGE: 
%  latentheat_melting = gsw_latentheat_melting(SA,p)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of melting.  It is defined in terms
%  of Absolute Salinity, SA, and sea pressure, p, and is valid in the  
%  ranges 0 < SA < 42 g kg^-1 and 0 < p < 10,000 dbar.  This is based on
%  the IAPWS Releases IAPWS-09 (for pure water), IAPWS-08 (for the saline
%  compoonent of seawater and IAPWS-06 for ice Ih.  
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & p need to have the same dimensions.
%
% OUTPUT:
%  latentheat_melting  =  latent heat of melting                   [ J/kg ]     
%
% AUTHOR:  
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
% IAPWS, 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic 
%  Properties of Seawater. The International Association for the Properties 
%  of Water and Steam. Berlin, Germany, September 2008.  This Release is 
%  known as IAPWS-09.  
%
% IAPWS, 2009a: Revised Release on the Equation of State 2006 for H2O Ice 
%  Ih. The International Association for the Properties of Water and Steam.
%  Doorwerth, The Netherlands, September 2009. This Release is known as 
%  IAPWS-06
%
% IAPWS, 2009b: Supplementary Release on a Computationally Efficient 
%  Thermodynamic Formulation for Liquid Water for Oceanographic Use. The 
%  International Association for the Properties of Water and Steam. 
%  Doorwerth, The Netherlands, September 2009.  This Release is known as 
%  IAPWS-09.  
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.34 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_latentheat_melting:  Requires two input arguments')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp ~= ms | np ~= ns)
    error('gsw_latentheat_melting: SA and p must have same dimensions')
end

if ms == 1
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

saturation_fraction = zeros(size(SA));

tf = gsw_t_freezing(SA,p,saturation_fraction);
latentheat_melting = 1000.*(gsw_chem_potential_water_t_exact(SA,tf,p) ... 
                     - (gsw_T0 + tf).*gsw_t_deriv_chem_potential_water_t_exact(SA,tf,p)) ... 
                     - gsw_enthalpy_ice(tf,p);

if transposed
    latentheat_melting = latentheat_melting.';
end

end
