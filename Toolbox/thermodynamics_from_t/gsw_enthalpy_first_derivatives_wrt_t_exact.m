function [h_SA_wrt_t, h_T_wrt_t, h_P_wrt_t] = gsw_enthalpy_first_derivatives_wrt_t_exact(SA,t,p)

% gsw_enthalpy_first_derivatives_wrt_t_exact           first derivatives of
%                                                                  enthalpy
%==========================================================================
%
% USAGE:
%  [h_SA_wrt_t, h_T_wrt_t, h_P_wrt_t] = gsw_enthalpy_first_derivatives_wrt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the following three derivatives of specific enthalpy, h.
%  These derivatives are done with respect to in-situ temperature t (in the
%  case of h_T_wrt_t) or at constant in-situ tempertature (in the cases of
%  h_SA_wrt_t and h_P_wrt_t).  
%   (1) h_SA_wrt_t, the derivative with respect to Absolute Salinity at 
%       constant t and p.
%   (2) h_T_wrt_t, derivative with respect to in-situ temperature t at 
%       constant SA and p. 
%   (3) h_P_wrt_t, derivative with respect to pressure P (in Pa) at constant 
%       SA and t.  This output has the same dimensions as specific volume,
%       but is not equal to specific volume.  
%
%  Note that this function uses the full Gibbs function.  This function
%  avoids the Nan that would exist in h_sub_SA at SA=0 if it were
%  evaluated in the straightforward way from the gibbs function (as in the
%  commented line 111 of the code below).  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  h_SA_wrt_t =  The first derivative of specific enthalpy with respect to 
%                Absolute Salinity at constant t and p.     
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_T_wrt_t  =  The first derivative of specific enthalpy with respect to 
%                in-situ temperature, t, at constant SA and p. [ J/(kg K) ]
%
%  h_P_wrt_t  =  The first derivative of specific enthalpy with respect to 
%                pressure P (in Pa) at constant SA and t.        [ m^3/kg ]
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
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_enthalpy_first_derivatives_wrt_t_exact:  Requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
    error('gsw_enthalpy_first_derivatives_wrt_t_exact: SA and t do not have the same dimensions')
end

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA.
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_enthalpy_first_derivatives_wrt_t_exact: The dimensions of p do not agree')
end

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0;

h_T_wrt_t = gsw_cp_t_exact(SA,t,p);
h_P_wrt_t = gsw_specvol_t_exact(SA,t,p) - (273.15 + t).*gsw_gibbs(0,1,1,SA,t,p);

%h_SA_wrt_t = gsw_gibbs(1,0,0,SA,t,p) - (273.15 + t).*gsw_gibbs(1,1,0,SA,t,p);

sfac = 0.0248826675584615; % sfac = 1/(40*(35.16504/35)).
x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*1e-4; %Note.The input pressure (p) is sea pressure in units of dbar.
    
    g08_SA = 8645.36753595126 + z.*(-6620.98308089678 + ...
        z.*(769.588305957198 + z.*(-193.0648640214916 + (31.6816345533648 - 5.24960313181984.*z).*z))) + ...
        x.*(-7296.43987145382 + x.*(8103.20462414788 + ...
        y.*(2175.341332000392 + y.*(-274.2290036817964 + ...
        y.*(197.4670779425016 + y.*(-68.5590309679152 + 9.98788038278032.*y))) - 90.6734234051316.*z) + ...
        x.*(-5458.34205214835 - 980.14153344888.*y + ...
        x.*(2247.60742726704 - 340.1237483177863.*x + 220.542973797483.*y) + 180.142097805543.*z) + ...
        z.*(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644.*z).*z)) + ...
        z.*(598.378809221703 + z.*(-156.8822727844005 + (204.1334828179377 - 10.23755797323846.*z).*z)) + ...
        y.*(-1480.222530425046 + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-129.1994027934126 + z.*(1149.174198007428 + z.*(-162.5751787551336 + 76.9195462169742.*z)) + ...
        y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(2.626801985426835 + 703.695562834065.*z))))) + ...
        y.*(1187.3715515697959 + z.*(1458.233059470092 + ...
        z.*(-687.913805923122 + z.*(249.375342232496 + z.*(-63.313928772146 + 14.09317606630898.*z)))) + ...
        y.*(1760.062705994408 + y.*(-450.535298526802 + ...
        y.*(182.8520895502518 + y.*(-43.3206481750622 + 4.26033941694366.*y) + ...
        z.*(-595.457483974374 + (149.452282277512 - 72.9745838003176.*z).*z)) + ...
        z.*(1388.489628266536 + z.*(-409.779283929806 + (227.123395681188 - 22.2565468652826.*z).*z))) + ...
        z.*(-1721.528607567954 + z.*(674.819060538734 + ...
        z.*(-356.629112415276 + (88.4080716616 - 15.84003094423364.*z).*z)))));

    gibbs_SA = 0.5.*sfac.*g08_SA;   
    
    g08_SA_T = 1187.3715515697959 + z.*(1458.233059470092 + ...
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
         
    gibbs_SA_T = 0.5.*sfac.*0.025.*g08_SA_T;

h_SA_wrt_t = gibbs_SA - (273.15 + t).*gibbs_SA_T;

if transposed
    h_T_wrt_t = h_T_wrt_t.';
    h_SA_wrt_t = h_SA_wrt_t.';
    h_P_wrt_t = h_P_wrt_t.';
end

end