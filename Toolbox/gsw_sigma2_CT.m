function sigma2 = gsw_sigma2_CT(SA,CT)

% gsw_sigma2                       potential density anomaly with reference
%                              sea pressure of 2000 dbar (48-term equation)
%==========================================================================
% This function has changed name to "gsw_sigma2"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sigma2_CT"
%==========================================================================

warning('This function has changed name to "gsw_sigma2".  This is the final version of GSW that will support "gsw_sigma2_CT" ')

sigma2 = gsw_sigma2(SA,CT);

end
