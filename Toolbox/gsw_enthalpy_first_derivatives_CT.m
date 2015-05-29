function [h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT(SA,CT,p)

% gsw_enthalpy_first_derivatives_CT           first derivatives of enthalpy
%==========================================================================
% This function has changed name to "gsw_enthalpy_first_derivatives"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_enthalpy_first_derivatives_CT"
%==========================================================================

warning('This function has changed name to "gsw_enthalpy_first_derivatives".  This is the final version of GSW that will support "gsw_enthalpy_first_derivatives_CT" ')

[h_SA, h_CT] = gsw_enthalpy_first_derivatives(SA,CT,p);


end