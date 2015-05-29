function sigma1 = gsw_sigma1_CT(SA,CT)

% gsw_sigma1                       potential density anomaly with reference
%                              sea pressure of 1000 dbar (48-term equation)
%==========================================================================
% This function has changed name to "gsw_sigma1"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sigma1_CT"
%==========================================================================

warning('This function has changed name to "gsw_sigma1".  This is the final version of GSW that will support "gsw_sigma1_CT" ')

sigma1 = gsw_sigma1(SA,CT);

end
