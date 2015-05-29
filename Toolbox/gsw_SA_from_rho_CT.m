function SA = gsw_SA_from_rho_CT(rho,CT,p)

% gsw_SA_from_rho                            Absolute Salinity from density
%==========================================================================
% This function has changed name to "gsw_SA_from_rho"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_SA_from_rho_CT"
%==========================================================================

warning('This function has changed name to "gsw_SA_from_rho".  This is the final version of GSW that will support "gsw_SA_from_rho_CT" ')

SA = gsw_SA_from_rho(rho,CT,p);

end
