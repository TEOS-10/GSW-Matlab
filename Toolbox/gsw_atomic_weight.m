function atomic_weight = gsw_atomic_weight

% gsw_atomic_weight                 mole-weighted atomic weight of sea salt
%==========================================================================
%
% USAGE:
%  atomic_weight = gsw_atomic_weight
%
% DESCRIPTION:
%  This function returns the mole-weighted atomic weight of sea salt of
%  Reference Composition, which is 31.4038218 g/mol.  This has been
%  defined as part of the Reference-Composition Salinity Scale of 2008 
%  (Millero et al., 2008). 
%
% OUTPUT:
%  atomic_weight = mole-weighted atomic weight of sea salt of Reference 
%                  Composition                                    [ g/mol ]
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
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%     See Eqn. (5.3) of this paper.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

atomic_weight = 31.4038218;

end
