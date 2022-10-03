function valence_factor = gsw_valence_factor

% gsw_valence_factor          valence factor of Reference Salinity sea salt
%==========================================================================
%
% USAGE:
%  valence_factor = gsw_valence_factor
%
% DESCRIPTION:
%  This function returns the valence factor of sea salt of Reference 
%  Composition, 1.2452898.  This valence factor is exact, and follows from
%  the definition of the Reference-Composition Salinity Scale 2008 of
%  Millero et al. (2008).  The valence factor is the mole-weighted square
%  of the charges, Z, of the ions comprising Reference Composition sea salt.  
%
% OUTPUT:
%  valence_factor  =  valence factor of sea salt of Reference Composition
%                                                              [ unitless ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Table D.4 of this TEOS-10 Manual. 
%
%  Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%     See Eqn. (5.9) of this paper. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

valence_factor = 1.2452898;

end
