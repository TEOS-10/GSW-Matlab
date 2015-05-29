function specvol_anom_t_exact = gsw_specvol_anom_t_exact(SA,t,p)

% gsw_specvol_anom_standard_t_exact                 specific volume anomaly
%==========================================================================
% This function has changed name to gsw_specvol_anom_standard_t_exact
% Note that this is the last version of GSW (Version 3.05) that will 
% contain the function "gsw_specvol_anom_t_exact".  Change your code to use
% "gsw_specvol_anom_standard_t_exact".
%==========================================================================

warning('This function has changed name to "gsw_specvol_anom_standard_t_exact". This is the last version of GSW that will support "gsw_specvol_anom_t_exact"')

specvol_anom_t_exact = gsw_specvol_anom_standard_t_exact(SA,t,p);

end
