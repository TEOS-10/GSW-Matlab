function thermobaric = gsw_thermobaric_CT(SA,CT,p)

% gsw_thermobaric                thermobaric coefficient (48-term equation)
%==========================================================================
% This function has changed name to "gsw_thermobaric"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_thermobaric_CT"
%==========================================================================

warning('This function has changed name to "gsw_thermobaric".  This is the final version of GSW that will support "gsw_thermobaric_CT" ')

thermobaric = gsw_thermobaric(SA,CT,p);

end
