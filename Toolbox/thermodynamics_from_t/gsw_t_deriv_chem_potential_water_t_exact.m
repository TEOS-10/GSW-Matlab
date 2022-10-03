function chem_potential_water_dt  = gsw_t_deriv_chem_potential_water_t_exact(SA,t,p)

% gsw_t_deriv_chem_potential_water_t_exact       the temperature derivative
%                                of chemical potential of water in seawater
%==========================================================================
%
% USAGE:
%  chem_potential_water_dt = gsw_t_deriv_chem_potential_water_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the temperature derivative of the chemical potential of water
%  in seawater so that it is valid at exactly SA = 0.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  chem_potential_water_dt  =  temperature derivative of the chemical 
%                           potential of water in seawater  [ J g^-1 K^-1 ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_t_deriv_chem_potential_water_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_t_deriv_chem_potential_water_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_t_deriv_chem_potential_water_t_exact: Inputs array dimensions arguments do not agree')
end %if

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

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
db2Pa = 1e-4;

x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*db2Pa; %Note that the input pressure (p) is sea pressure in units of dbar.

g03_T = 5.90578347909402 + z.*(-270.983805184062 + ...
        z.*(776.153611613101 + z.*(-196.51255088122 + (28.9796526294175 - 2.13290083518327.*z).*z))) + ...
        y.*(-24715.571866078 + z.*(2910.0729080936 + ...
        z.*(-1513.116771538718 + z.*(546.959324647056 + z.*(-111.1208127634436 + 8.68841343834394.*z)))) + ...
        y.*(2210.2236124548363 + z.*(-2017.52334943521 + ...
        z.*(1498.081172457456 + z.*(-718.6359919632359 + (146.4037555781616 - 4.9892131862671505.*z).*z))) + ...
        y.*(-592.743745734632 + z.*(1591.873781627888 + ...
        z.*(-1207.261522487504 + (608.785486935364 - 105.4993508931208.*z).*z)) + ...
        y.*(290.12956292128547 + z.*(-973.091553087975 + ...
        z.*(602.603274510125 + z.*(-276.361526170076 + 32.40953340386105.*z))) + ...
        y.*(-113.90630790850321 + y.*(21.35571525415769 - 67.41756835751434.*z) + ...
        z.*(381.06836198507096 + z.*(-133.7383902842754 + 49.023632509086724.*z)))))));
    
g08_T = x2.*(168.072408311545 + ...
        x.*(-493.407510141682 + x.*(543.835333000098 + x.*(-196.028306689776 + 36.7571622995805.*x) + ...
        y.*(-137.1145018408982 + y.*(148.10030845687618 + y.*(-68.5590309679152 + 12.4848504784754.*y))) - ...
        22.6683558512829.*z) + z.*(-175.292041186547 + (83.1923927801819 - 29.483064349429.*z).*z) + ...
        y.*(-86.1329351956084 + z.*(766.116132004952 + z.*(-108.3834525034224 + 51.2796974779828.*z)) + ...
        y.*(-30.0682112585625 - 1380.9597954037708.*z + y.*(3.50240264723578 + 938.26075044542.*z)))));
    
g08_SA_T = 1187.3715515697959 + ...
        x.*(-1480.222530425046 + x.*(2175.341332000392 + x.*(-980.14153344888 + 220.542973797483.*x) + ...
        y.*(-548.4580073635929 + y.*(592.4012338275047 + y.*(-274.2361238716608 + 49.9394019139016.*y))) - ...
        90.6734234051316.*z) + z.*(-525.876123559641 + (249.57717834054571 - 88.449193048287.*z).*z) + ...
        y.*(-258.3988055868252 + z.*(2298.348396014856 + z.*(-325.1503575102672 + 153.8390924339484.*z)) + ...
        y.*(-90.2046337756875 - 4142.8793862113125.*z + y.*(10.50720794170734 + 2814.78225133626.*z))));
    
kg2g = 1e-3;
chem_potential_water_dt =  kg2g*((g03_T + g08_T).*0.025 - 0.5.*sfac.*0.025.*SA.*g08_SA_T);

% Note. The kg2g, a factor of 1e-3, is needed to convert the output of this function into units of J/g. 
% See section (2.9) of the TEOS-10 Manual.

if transposed
    chem_potential_water_dt = chem_potential_water_dt.';
end

end
