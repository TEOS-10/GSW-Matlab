function rho_CT = gsw_rho_CT(SA,CT,p)

% gsw_rho_CT                             in-situ density (48-term equation)
%==========================================================================
% This function has changed name to "gsw_rho"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_rho_CT"
%==========================================================================

warning('This function has changed name to "gsw_rho".  This is the final version of GSW that will support "gsw_rho_CT" ')

rho_CT = gsw_rho(SA,CT,p);

end

