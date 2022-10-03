function [pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA,CT)

% gsw_pt_second_derivatives     second derivatives of potential temperature 
% =========================================================================
%
% USAGE:
%  [pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of potential 
%  temperature (the regular potential temperature which has a reference 
%  sea pressure of 0 dbar), 
%   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity 
%       at constant Conservative Temperature,
%   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
%       and Absolute Salinity, and
%   (3) pt_CT_CT, the second derivative with respect to Conservative 
%       Temperature at constant Absolute Salinity. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt_SA_SA  =  The second derivative of potential temperature (the 
%               regular potential temperature which has reference sea 
%               pressure of 0 dbar) with respect to Absolute Salinity 
%               at constant Conservative Temperature.  
%               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  pt_SA_CT  =  The derivative of potential temperature with respect 
%               to Absolute Salinity and Conservative Temperature.   
%               pt_SA_CT has units of:                         [ 1/(g/kg) ]
%  pt_CT_CT  =  The second derivative of potential temperature (the 
%               regular one with p_ref = 0 dbar) with respect to 
%               Conservative Temperature at constant SA.  
%               pt_CT_CT has units of:                              [ 1/K ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (7th September, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org. 
%    See Eqns. (A.12.9) and (A.12.10) of this TEOS-10 Manual.     
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_pt_second_derivatives:  Requires two inputs')
end %if

if ~(nargout == 3)
   error('gsw_pt_second_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_pt_second_derivatives:  SA and CT must have same dimensions')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

dSA = 1e-3;                 % increment of Absolute Salinity is 0.001 g/kg.

SA_l = nan(size(SA));
SA_l(SA>=dSA) = SA(SA>=dSA) - dSA;
SA_l(SA<dSA) = 0;

SA_u = SA + dSA;

[pt_SA_l, pt_CT_l] = gsw_pt_first_derivatives(SA_l,CT);
[pt_SA_u, pt_CT_u] = gsw_pt_first_derivatives(SA_u,CT);

pt_SA_SA = (pt_SA_u - pt_SA_l)./(SA_u - SA_l);
% pt_SA_CT = (pt_CT_u - pt_CT_l)./(SA_u - SA_l); % can calculate this either way;

dCT  = 1e-2;     % increment of Conservative Temperature is 0.01 degrees C;
rec_2dCT = 50;   % rec_2dCT = 1./(2.*dCT);
CT_l = CT - dCT;
CT_u = CT + dCT;

[pt_SA_l, pt_CT_l] = gsw_pt_first_derivatives(SA,CT_l);
[pt_SA_u, pt_CT_u] = gsw_pt_first_derivatives(SA,CT_u);

pt_SA_CT = (pt_SA_u - pt_SA_l).*rec_2dCT;
pt_CT_CT = (pt_CT_u - pt_CT_l).*rec_2dCT;

if transposed   
    pt_SA_CT = pt_SA_CT.';
    pt_CT_CT = pt_CT_CT.';
end

end
