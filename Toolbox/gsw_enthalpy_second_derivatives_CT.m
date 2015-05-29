function [h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives_CT(SA,CT,p)

% gsw_enthalpy_second_derivatives_CT         second derivatives of enthalpy
%                                                        (48-term equation)
%==========================================================================
% This function has changed name to "gsw_enthalpy_second_derivatives"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_enthalpy_second_derivatives_CT"
%==========================================================================

warning('This function has changed name to "gsw_enthalpy_second_derivatives".  This is the final version of GSW that will support "gsw_enthalpy_second_derivatives_CT" ')

[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA,CT,p);

end
