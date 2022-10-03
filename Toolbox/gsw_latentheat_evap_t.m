function latentheat_evap = gsw_latentheat_evap_t(SA,t)

% gsw_latentheat_evap_t                          latent heat of evaporation 
%==========================================================================
%
% USAGE: 
%  latentheat_evap = gsw_latentheat_evap_t(SA,t)
%
% DESCRIPTION:
%  Calculates latent heat, or enthalpy, of evaporation at p = 0 (the 
%  surface).  It is defined as a function of Absolute Salinity, SA, and
%  in-situ temperature, t, and is valid in the ranges 0 < SA < 40 g/kg 
%  and 0 < CT < 42 deg C. The errors range between -0.4 and 0.6 J/kg.
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  latentheat_evap = latent heat of evaporation                    [ J/kg ]
%
% AUTHOR:  
%  Paul Barker, Trevor McDougall & Rainer Feistel      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See section 3.39 of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_latentheat_evap_t: Requires two input arguments')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);

if (mt ~= ms | nt ~= ns)
    error('gsw_latentheat_evap_t: SA and t must have same dimensions')
end

if ms == 1
    SA = SA.';
    t = t.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
 
CT = gsw_CT_from_pt(SA,t);

latentheat_evap = gsw_latentheat_evap_CT(SA,CT);

if transposed
    latentheat_evap = latentheat_evap.';
end

end
