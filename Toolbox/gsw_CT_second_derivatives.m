function [CT_SA_SA, CT_SA_pt, CT_pt_pt] = gsw_CT_second_derivatives(SA,pt)

% gsw_CT_second_derivatives                           second derivatives of  
%                                                  Conservative Temperature 
%==========================================================================
%
% USAGE:
%  [CT_SA_SA, CT_SA_pt, CT_pt_pt] = gsw_CT_second_derivatives(SA,pt)
%
% DESCRIPTION:
%  Calculates the following three, second-order derivatives of Conservative 
%  Temperature
%   (1) CT_SA_SA, the second derivative with respect to Absolute Salinity  
%       at constant potential temperature (with p_ref = 0 dbar),
%   (2) CT_SA_pt, the derivative with respect to potential temperature
%       (the regular potential temperature which is referenced to 0 dbar)
%       and Absolute Salinity, and
%   (3) CT_pt_pt, the second derivative with respect to potential 
%       temperature (the regular potential temperature which is referenced 
%       to 0 dbar) at constant Absolute Salinity. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]   
%         (whose reference pressure is 0 dbar)
%
%  SA & pt need to have the same dimensions.
%
% OUTPUT:
%  CT_SA_SA  =  The second derivative of Conservative Temperature with 
%               respect to Absolute Salinity at constant potential 
%               temperature (the regular potential temperature which 
%               has reference sea pressure of 0 dbar).  
%               CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  CT_SA_pt  =  The derivative of Conservative Temperature with 
%               respect to potential temperature (the regular one with 
%               p_ref = 0 dbar) and Absolute Salinity.   
%               CT_SA_pt has units of:                        [ 1/(g/kg) ]
%  CT_pt_pt  =  The second derivative of Conservative Temperature with 
%               respect to potential temperature (the regular one with 
%               p_ref = 0 dbar) at constant SA.   
%               CT_pt_pt has units of:                              [ 1/K ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (7th September, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See appendix A.12 of this TEOS-10 Manual.    
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_CT_second_derivatives:  Requires two inputs')
end

if ~(nargout == 3)
   error('gsw_CT_second_derivatives:  Requires three outputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_CT_second_derivatives:  SA and pt must have same dimensions')
end

if ms == 1
    SA = SA.';
    pt = pt.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

dSA = 1e-3;                 % increment of Absolute Salinity is 0.001 g/kg.
SA_l = nan(size(SA));
SA_l(SA >= dSA) = SA(SA >= dSA) - dSA;
SA_l(SA < dSA) = 0;

SA_u = SA + dSA;

[CT_SA_l, dummy] = gsw_CT_first_derivatives(SA_l,pt);
[CT_SA_u, dummy] = gsw_CT_first_derivatives(SA_u,pt);

CT_SA_SA = nan(size(SA));
CT_SA_SA(SA_u ~= SA_l) = (CT_SA_u(SA_u ~= SA_l) - CT_SA_l(SA_u ~= SA_l))./ ...
                         (SA_u(SA_u ~= SA_l) - SA_l(SA_u ~= SA_l));

dpt  = 1e-2;        % increment of potential temperature is 0.01 degrees C.
rec_2dpt = 50;      % rec_2dpt = 1./(2.*dpt);
pt_l = pt - dpt;
pt_u = pt + dpt;

[CT_SA_l, CT_pt_l] = gsw_CT_first_derivatives(SA,pt_l);
[CT_SA_u, CT_pt_u] = gsw_CT_first_derivatives(SA,pt_u);

CT_SA_pt = (CT_SA_u - CT_SA_l).*rec_2dpt;
CT_pt_pt = (CT_pt_u - CT_pt_l).*rec_2dpt;

if transposed
    CT_SA_SA = CT_SA_SA.';
    CT_SA_pt = CT_SA_pt.';
    CT_pt_pt = CT_pt_pt.';
end

end
