function sigma3 = gsw_sigma3_CT(SA,CT)

% gsw_sigma3                       potential density anomaly with reference
%                              sea pressure of 3000 dbar (48-term equation)
%==========================================================================
% This function has changed name to "gsw_sigma3"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sigma3_CT"
%==========================================================================

warning('This function has changed name to "gsw_sigma3".  This is the final version of GSW that will support "gsw_sigma3_CT" ')

sigma3 = gsw_sigma3(SA,CT);

end
