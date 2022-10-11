function adiabatic_lapse_rate = gsw_adiabatic_lapse_rate_t_exact(SA,t,p)

% gsw_adiabatic_lapse_rate_from_t                      adiabatic lapse rate
%==========================================================================
%
% This function has changed its name, it is now called 
% gsw_adiabatic_lapse_rate_from_t.
%
% USAGE:
%  adiabatic_lapse_rate = gsw_adiabatic_lapse_rate_from_t(SA,t,p)
%
% DESCRIPTION:
%  Calculates the adiabatic lapse rate of sea water
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  adiabatic_lapse_rate  =  adiabatic lapse rate                   [ K/Pa ]
%    Note.  The output is in unit of degress Celsius per Pa,
%      (or equivilently K/Pa) not in units of K/dbar. 
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqn. (2.22.1) of this TEOS-10 Manual.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

warning('The function "gsw_adiabatic_lapse_rate_t_exact" has changed to "gsw_adiabatic_lapse_rate_from_t"');

adiabatic_lapse_rate = gsw_adiabatic_lapse_rate_from_t(SA,t,p);

end
