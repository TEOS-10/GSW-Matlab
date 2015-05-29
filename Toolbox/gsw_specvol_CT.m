function specvol = gsw_specvol_CT(SA,CT,p)

% gsw_specvol                            specific volume (48-term equation)
%==========================================================================
% This function has changed name to "gsw_specvol"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_specvol_CT"
%==========================================================================

warning('This function has changed name to "gsw_specvol".  This is the final version of GSW that will support "gsw_specvol_CT" ')

specvol = gsw_specvol(SA,CT,p);

end
