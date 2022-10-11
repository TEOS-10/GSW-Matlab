function chem_potential_water_t_exact = gsw_chem_potential_water_t_exact(SA,t,p)

% gsw_chem_potential_water_t_exact              chemical potential of water
%                                                               in seawater
%==========================================================================
%
% USAGE:
%  chem_potential_water_t_exact  = gsw_chem_potential_water_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the chemical potential of water in seawater.
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
%  chem_potential_water_t_exact  =  chemical potential of water in seawater
%                                                                   [ J/g ]
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
   error('gsw_chem_potential_water_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_chem_potential_water_t_exact: SA and t must have same dimensions')
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
    error('gsw_chem_potential_water_t_exact: Inputs array dimensions arguments do not agree')
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

x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*1e-4; %Note that the input pressure (p) is sea pressure in units of dbar.

g03_g = 101.342743139674 + z.*(100015.695367145 + ...
    z.*(-2544.5765420363 + z.*(284.517778446287 + ...
    z.*(-33.3146754253611 + (4.20263108803084 - 0.546428511471039.*z).*z)))) + ...
    y.*(5.90578347909402 + z.*(-270.983805184062 + ...
    z.*(776.153611613101 + z.*(-196.51255088122 + (28.9796526294175 - 2.13290083518327.*z).*z))) + ...
    y.*(-12357.785933039 + z.*(1455.0364540468 + ...
    z.*(-756.558385769359 + z.*(273.479662323528 + z.*(-55.5604063817218 + 4.34420671917197.*z)))) + ...
    y.*(736.741204151612 + z.*(-672.50778314507 + ...
    z.*(499.360390819152 + z.*(-239.545330654412 + (48.8012518593872 - 1.66307106208905.*z).*z))) + ...
    y.*(-148.185936433658 + z.*(397.968445406972 + ...
    z.*(-301.815380621876 + (152.196371733841 - 26.3748377232802.*z).*z)) + ...
    y.*(58.0259125842571 + z.*(-194.618310617595 + ...
    z.*(120.520654902025 + z.*(-55.2723052340152 + 6.48190668077221.*z))) + ...
    y.*(-18.9843846514172 + y.*(3.05081646487967 - 9.63108119393062.*z) + ...
    z.*(63.5113936641785 + z.*(-22.2897317140459 + 8.17060541818112.*z))))))));

g08_g = x2.*(1416.27648484197 + ...
    x.*(-2432.14662381794 + x.*(2025.80115603697 + ...
    y.*(543.835333000098 + y.*(-68.5572509204491 + ...
    y.*(49.3667694856254 + y.*(-17.1397577419788 + 2.49697009569508.*y))) - 22.6683558512829.*z) + ...
    x.*(-1091.66841042967 - 196.028306689776.*y + ...
    x.*(374.60123787784 - 48.5891069025409.*x + 36.7571622995805.*y) + 36.0284195611086.*z) + ...
    z.*(-54.7919133532887 + (-4.08193978912261 - 30.1755111971161.*z).*z)) + ...
    z.*(199.459603073901 + z.*(-52.2940909281335 + (68.0444942726459 - 3.41251932441282.*z).*z)) + ...
    y.*(-493.407510141682 + z.*(-175.292041186547 + (83.1923927801819 - 29.483064349429.*z).*z) + ...
    y.*(-43.0664675978042 + z.*(383.058066002476 + z.*(-54.1917262517112 + 25.6398487389914.*z)) + ...
    y.*(-10.0227370861875 - 460.319931801257.*z + y.*(0.875600661808945 + 234.565187611355.*z))))) + ...
    y.*(168.072408311545));

g_SA_part = 8645.36753595126 + ...
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
    y.*(1187.3715515697959);

g2kg = 1e-3;

chem_potential_water_t_exact =  g2kg*(g03_g + g08_g - 0.5.*sfac.*SA.*g_SA_part);
% Note. The g2kg, a factor of 1e-3, is needed to convert the output of this function into units of J/g. 
% See section (2.9) of the TEOS-10 Manual.

if transposed
    chem_potential_water_t_exact = chem_potential_water_t_exact.';
end

end
