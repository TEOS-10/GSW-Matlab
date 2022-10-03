function pot_enthalpy_ice = gsw_pot_enthalpy_from_specvol_ice(specvol_ice,p)

% gsw_pot_enthalpy_from_specvol_ice                 potential enthalpy from
%                                                    specific volume of ice                              
%==========================================================================
%
% USAGE:
%  pot_enthalpy_ice = gsw_pot_enthalpy_from_specvol_ice(specvol_ice,p)
%
% DESCRIPTION:
%  Calculates the potential enthalpy of ice from the specific volume
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar. 
%
% INPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OUTPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
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
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_pot_enthalpy_from_specvol_ice:  Requires two inputs')
end

[mh,nh] = size(specvol_ice);
[mp,np] = size(p);

if (mp == 1) & (np == 1)           
    p = p*ones(mh,nh);
elseif (nh == np) & (mp == 1)         
    p = p(ones(1,mh), :);             
elseif (mh == mp) & (np == 1)        
    p = p(:,ones(1,nh));               
elseif (nh == mp) & (np == 1)        
    p = p.';                     
    p = p(ones(1,mh), :);          
elseif (mh == mp) & (nh == np)
    % ok
else
    error('gsw_pot_enthalpy_from_specvol_ice: Inputs array dimensions arguments do not agree')
end

if mh == 1
    specvol_ice = specvol_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

rho_ice = 1./specvol_ice;

t_ice = gsw_t_from_rho_ice(rho_ice,p);

pt0_ice = gsw_pt0_from_t_ice(t_ice,p);

pot_enthalpy_ice = gsw_pot_enthalpy_from_pt_ice(pt0_ice);

if transposed
   pot_enthalpy_ice = pot_enthalpy_ice.';
end

end
