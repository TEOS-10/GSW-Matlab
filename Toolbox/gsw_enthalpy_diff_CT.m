function enthalpy_diff = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff                 difference of enthalpy at two pressures
%                                                        (48-term equation)
%==========================================================================
% This function has changed name to "gsw_enthalpy_diff"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_enthalpy_diff_CT"
%==========================================================================

warning('This function has changed name to "gsw_enthalpy_diff".  This is the final version of GSW that will support "gsw_enthalpy_diff_CT" ')

enthalpy_diff = gsw_enthalpy_diff(SA,CT,p_shallow,p_deep);

end
