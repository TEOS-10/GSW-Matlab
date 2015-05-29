function alpha = gsw_alpha_CT(SA,CT,p)

% gsw_alpha                   thermal expansion coefficient with respect to 
%                               Conservative Temperature (48-term equation)
%==========================================================================
% This function has changed name to gsw_alpha
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_alpha_CT"
%==========================================================================

warning('This function has changed name to "gsw_alpha". This is the final version of GSW that will support "gsw_alpha_CT"')

alpha = gsw_alpha(SA,CT,p);

end
