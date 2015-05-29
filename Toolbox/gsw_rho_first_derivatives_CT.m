function [drho_dSA, drho_dCT, drho_dP] = gsw_rho_first_derivatives_CT(SA,CT,p)

% gsw_rho_first_derivatives                SA, CT and p partial derivatives
%                                             of density (48-term equation)
%==========================================================================
% This function has changed name to "gsw_rho_first_derivatives"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_rho_first_derivatives_CT"
%==========================================================================

warning('This function has changed name to "gsw_rho_first_derivatives".  This is the final version of GSW that will support "gsw_rho_first_derivatives_CT" ')

[drho_dSA, drho_dCT, drho_dP] = gsw_rho_first_derivatives(SA,CT,p);

end
