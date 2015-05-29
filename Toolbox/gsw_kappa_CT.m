function kappa = gsw_kappa_CT(SA,CT,p)

% gsw_kappa                   isentropic compressibility (48-term equation)     
%==========================================================================
% This function has changed name to "gsw_kappa"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_kappa_CT"
%==========================================================================

warning('This function has changed name to "gsw_kappa".  This is the final version of GSW that will support "gsw_kappa_CT" ')

kappa = gsw_kappa(SA,CT,p);

end