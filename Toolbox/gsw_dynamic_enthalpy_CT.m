function dynamic_enthalpy = gsw_dynamic_enthalpy_CT(SA,CT,p)

% gsw_dynamic_enthalpy                         dynamic enthalpy of seawater
%                                                        (48-term equation)
%==========================================================================
% This function has changed name to "gsw_dynamic_enthalpy"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_dynamic_enthalpy_CT"
%==========================================================================

warning('This function has changed name to "gsw_dynamic_enthalpy".  This is the final version of GSW that will support "gsw_dynamic_enthalpy_CT" ')

dynamic_enthalpy = gsw_dynamic_enthalpy(SA,CT,p);

end
