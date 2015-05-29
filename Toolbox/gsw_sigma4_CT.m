function sigma4 = gsw_sigma4_CT(SA,CT)

% gsw_sigma4                       potential density anomaly with reference
%                              sea pressure of 4000 dbar (48-term equation)
%==========================================================================
% This function has changed name to "gsw_sigma4"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sigma4_CT"
%==========================================================================

warning('This function has changed name to "gsw_sigma4".  This is the final version of GSW that will support "gsw_sigma4_CT" ')

sigma4 = gsw_sigma4(SA,CT);

end
