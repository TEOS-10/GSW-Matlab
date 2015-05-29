function [rho, alpha, beta] = gsw_rho_alpha_beta_CT(SA,CT,p)

% gsw_rho_alpha_beta            in-situ density, thermal expansion & saline 
%                                contraction coefficient (48-term equation)
%==========================================================================
% This function has changed name to "gsw_rho_alpha_beta"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_rho_alpha_beta_CT"
%==========================================================================

warning('This function has changed name to "gsw_rho_alpha_beta". This is the final version of GSW that will support "gsw_rho_alpha_beta"')

[rho, alpha, beta] = gsw_rho_alpha_beta(SA,CT,p);

end
