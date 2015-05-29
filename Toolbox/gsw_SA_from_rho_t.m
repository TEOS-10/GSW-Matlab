function SA = gsw_SA_from_rho_t(rho,t,p)

% gsw_SA_from_rho_t                          Absolute Salinity from density
% =========================================================================
%
% USAGE:
%  SA = gsw_SA_from_rho_t(rho,t,p)
% 
% DESCRIPTION:
%  Calculates the Absolute Salinity of a seawater sample, for given values
%  of its density, in-situ temperature and sea pressure (in dbar). 
%  This function uses the computationally-efficient 48-term expression for 
%  density in terms of SA, CT and p (McDougall et al., 2013).
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2013).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it. 
%     That is, it is 'density', not 'density anomaly'.
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & t are MxN.
%
% OUTPUT:
%  SA  =  Absolute Salinity.                                       [ g/kg ]
%   Note. This is expressed on the Reference-Composition Salinity
%     Scale of Millero et al. (2008). 
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.02 (15th November, 2012)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2013:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to J. Atm. Ocean. Technol., xx, yyy-zzz.
%
%  Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
   error('gsw_SA_from_rho_t:  Requires three inputs')
end %if

[md,nd] = size(rho);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= md | nt ~= nd)
    error('gsw_SA_from_rho_t: rho and t must have same dimensions')
end

if (mp == 1) & (np == 1)               % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)          % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)          % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                            % transposed then
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_SA_from_rho_t: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

v_lab = ones(size(rho))./rho;
CT_0 = gsw_CT_from_t(zeros(size(rho)),t,p);
v_0 = gsw_specvol(zeros(size(rho)),CT_0,p);
CT_50 = gsw_CT_from_t(50*ones(size(rho)),t,p);
v_50 = gsw_specvol(50*ones(size(rho)),CT_50,p);
 
SA = 50*(v_lab - v_0)./(v_50 - v_0);            % initial estimate of SA.
SA(SA < 0 | SA > 50) = NaN;

v_SA_t = 0.02.*(v_50 - v_0); %  initial estimate of the derivative of 
%                               specific volume with respect to
%                               Absolute Salinity at constant in-situ 
%                               temperature and pressure.  
                  
%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure 
%--------------------------------------------------------------------------

for Number_of_iterations = 1:2 
    SA_old = SA;
    CT_old = gsw_CT_from_t(SA_old,t,p);
    delta_v = gsw_specvol(SA_old,CT_old,p) - v_lab;
    SA = SA_old - delta_v./v_SA_t; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)  
    SA(SA < 0) = 0;
    SA_mean = 0.5*(SA + SA_old);
    CT_mean = gsw_CT_from_t(SA_mean,t,p);
    [rho,alpha,beta] = gsw_rho_alpha_beta(SA_mean,CT_mean,p);
    specvol = 1./rho;
    v_SA = - beta.*specvol; 
    v_CT = alpha.*specvol;
    dCT_dSA = gsw_dCT_dSA(SA_mean,t,p);% This is a call to get the 
%                                        derivative of CT with respect to
%                                        Absolute Salinity at constant 
%                                        in-situ temperature and pressure.  

    v_SA_t = v_SA + v_CT.*dCT_dSA;
%     v_SA_t is the derivative of specific volume with respect to
%     Absolute Salinity at constant in-situ temperature and pressure, 
%     while v_SA is the derivative of specific volume with respect to
%     Absolute Salinity at constant Conservative Temperature and pressure. 

    SA = SA_old - delta_v./v_SA_t;
    SA(SA < 0 | SA > 50) = NaN; 
end

% After two iterations of this modified Newton-Raphson iteration,
% the error in SA is no larger than 3x10^-12 g/kg, which 
% is machine precision for this calculation. 
 
if transposed
    SA = SA.';
end

end

%##########################################################################

function dCT_dSA = gsw_dCT_dSA(SA,t,p)
%--------------------------------------------------------------------------
% This function calculates the derivative of Conservative Temperature 
% with respect to Absolute Salinity at constant in-situ temperature and 
% pressure.  The equation for this derivative, in terms of the derivatives 
% of the Gibbs function is given in Eqn. (A.15.8) of the TEOS-10 Manual 
% (IOC et al., 2010).  The code below is identical to the first part of the
% GSW function gsw_beta_const_CT_t_exact. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%--------------------------------------------------------------------------

db2Pa = 1e-4;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
cp0 = 3991.86795711963;

pt0 = gsw_pt0_from_t(SA,t,p);

x2 = sfac.*SA;
x = sqrt(x2);
y = 0.025*t;
y_pt = 0.025*pt0;
z = db2Pa*p; %Note.The input pressure (p) is sea pressure in units of dbar.

g_SA_T_mod = 1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y.*(-548.4580073635929 + y.*(592.4012338275047 + y.*(-274.2361238716608 + 49.9394019139016.*y))) - ...
        90.6734234051316.*z) + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-258.3988055868252 + z.*(2298.348396014856 + z.*(-325.1503575102672 + 153.8390924339484.*z)) + ...
        y.*(-90.2046337756875 - 4142.8793862113125.*z + y.*(10.50720794170734 + 2814.78225133626.*z)))) + ...
        y.*(3520.125411988816 + y.*(-1351.605895580406 + ...
        y.*(731.4083582010072 + y.*(-216.60324087531103 + 25.56203650166196.*y) + ...
        z.*(-2381.829935897496 + (597.809129110048 - 291.8983352012704.*z).*z)) + ...
        z.*(4165.4688847996085 + z.*(-1229.337851789418 + (681.370187043564 - 66.7696405958478.*z).*z))) + ...
        z.*(-3443.057215135908 + z.*(1349.638121077468 + ...
        z.*(-713.258224830552 + (176.8161433232 - 31.68006188846728.*z).*z))));
g_SA_T_mod = 0.5*sfac*0.025*g_SA_T_mod;
   
g_SA_mod = 8645.36753595126 + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y_pt.*(2175.341332000392 + y_pt.*(-274.2290036817964 + ...
        y_pt.*(197.4670779425016 + y_pt.*(-68.5590309679152 + 9.98788038278032.*y_pt)))) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y_pt + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y_pt))) + ...
        y_pt.*(-1480.222530425046 +  y_pt.*(-129.1994027934126 + ...
        y_pt.*(-30.0682112585625 + y_pt.*(2.626801985426835 ))))) + ...
        y_pt.*(1187.3715515697959 + y_pt.*(1760.062705994408 + y_pt.*(-450.535298526802 + ...
        y_pt.*(182.8520895502518 + y_pt.*(-43.3206481750622 + 4.26033941694366.*y_pt)))));
g_SA_mod = 0.5*sfac*g_SA_mod;   

dCT_dSA = (g_SA_mod - (273.15 + pt0).*g_SA_T_mod)./cp0;

end
