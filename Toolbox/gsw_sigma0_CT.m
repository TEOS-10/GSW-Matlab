function sigma0 = gsw_sigma0_CT(SA,CT)

% gsw_sigma0                       potential density anomaly with reference
%                                 sea pressure of 0 dbar (48-term equation)
%==========================================================================
% This function has changed name to "gsw_sigma0"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sigma0_CT"
%==========================================================================

warning('This function has changed name to "gsw_sigma0".  This is the final version of GSW that will support "gsw_sigma0_CT" ')

sigma0 = gsw_sigma0(SA,CT);

end
