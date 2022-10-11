function [CT_SA, CT_pt] = gsw_CT_first_derivatives(SA,pt)

% gsw_CT_first_derivatives    first derivatives of Conservative Temperature
%==========================================================================
%
% USAGE:
%  [CT_SA, CT_pt] = gsw_CT_first_derivatives(SA,pt)
%
% DESCRIPTION:
%  Calculates the following two derivatives of Conservative Temperature
%  (1) CT_SA, the derivative with respect to Absolute Salinity at 
%      constant potential temperature (with pr = 0 dbar), and
%   2) CT_pt, the derivative with respect to potential temperature
%      (the regular potential temperature which is referenced to 0 dbar)
%      at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]   
%         (whose reference pressure is 0 dbar)
%
%  SA & pt need to have the same dimensions.
%
% OUTPUT:
%  CT_SA  =  The derivative of Conservative Temperature with respect to 
%            Absolute Salinity at constant potential temperature 
%            (the regular potential temperature which has reference 
%            sea pressure of 0 dbar).    
%            The CT_SA output has units of:                     [ K/(g/kg)]
%  CT_pt  =  The derivative of Conservative Temperature with respect to 
%            potential temperature (the regular one with pr = 0 dbar)
%            at constant SA. CT_pt is dimensionless.           [ unitless ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqns. (A.12.3a,b) and (A.15.8) of this TEOS-10 Manual.   
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_CT_first_derivatives:  Requires two inputs')
end %if

if ~(nargout == 2)
   error('gsw_CT_first_derivatives:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_CT_first_derivatives: SA and t must have same dimensions')
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

cp0 = gsw_cp0;             % from Eqn. 3.3.3 of IOC et al. (2010).
abs_pt = gsw_T0 + pt; 

CT_pt = - (abs_pt.*gsw_gibbs_pt0_pt0(SA,pt))./cp0;

%--------------------------------------------------------------------------

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
x2 = sfac.*SA;
x = sqrt(x2);
y_pt = 0.025*pt;

g_SA_T_mod = 1187.3715515697959 + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y_pt.*(-548.4580073635929 + y_pt.*(592.4012338275047 + y_pt.*(-274.2361238716608 + 49.9394019139016.*y_pt)))) + ...
        y_pt.*(-258.3988055868252 + y_pt.*(-90.2046337756875 + y_pt.*10.50720794170734))) + ...
        y_pt.*(3520.125411988816  + y_pt.*(-1351.605895580406 + ...
        y_pt.*(731.4083582010072  + y_pt.*(-216.60324087531103 + 25.56203650166196.*y_pt))));
g_SA_T_mod = 0.5*sfac*0.025*g_SA_T_mod;
   
g_SA_mod = 8645.36753595126 + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y_pt.*(2175.341332000392 + y_pt.*(-274.2290036817964 + ...
        y_pt.*(197.4670779425016 + y_pt.*(-68.5590309679152 + 9.98788038278032.*y_pt)))) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y_pt + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y_pt))) + ...
        y_pt.*(-1480.222530425046 + ...
        y_pt.*(-129.1994027934126 + ...
        y_pt.*(-30.0682112585625 + y_pt.*(2.626801985426835 ))))) + ...
        y_pt.*(1187.3715515697959 + ...
        y_pt.*(1760.062705994408 + y_pt.*(-450.535298526802 + ...
        y_pt.*(182.8520895502518 + y_pt.*(-43.3206481750622 + 4.26033941694366.*y_pt)))));
g_SA_mod = 0.5*sfac*g_SA_mod;   

CT_SA = (g_SA_mod - abs_pt.*g_SA_T_mod)./cp0;

if transposed
    CT_SA = CT_SA.';
    CT_pt = CT_pt.';
end

end
