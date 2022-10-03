function dynamic_enthalpy = gsw_dynamic_enthalpy(SA,CT,p)

% gsw_dynamic_enthalpy                         dynamic enthalpy of seawater
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  dynamic_enthalpy = gsw_dynamic_enthalpy(SA,CT,p)
%
% DESCRIPTION:
%  Calculates dynamic enthalpy of seawater using the computationally-
%  efficient expression for specific volume in terms of SA, CT and p
%  (Roquet et al., 2015).  Dynamic enthalpy is defined as enthalpy minus
%  potential enthalpy (Young, 2010). 
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  dynamic_enthalpy  =  dynamic enthalpy                           [ J/kg ]
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
%    See section 3.2 of this TEOS-10 Manual. 
%
%  McDougall, T. J., 2003: Potential enthalpy: A conservative oceanic 
%   variable for evaluating heat content and heat fluxes. Journal of 
%   Physical Oceanography, 33, 945-963.  
%    See Eqns. (18) and (22)
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  Young, W.R., 2010: Dynamic enthalpy, Conservative Temperature, and the
%   seawater Boussinesq approximation. Journal of Physical Oceanography, 
%   40, 394-400.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_dynamic_enthalpy: requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_dynamic_enthalpy: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
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
    error('gsw_dynamic_enthalpy: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
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

% % Set lower temperature limit for water that is much colder than the
% % freezing termperature.
% CT_frozen = - 1.79e-02 - 6.7157e-2.*SA - 9.2708e-04.*p;
% 
% if any(CT(:) < (CT_frozen(:) - 0.2))
%     [Icold] = find(CT < CT_frozen - 0.2);
%     CT(Icold) = NaN;
% end

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z = p.*1e-4;

h001 =  1.0769995862e-3; 
h002 = -3.0399571905e-5; 
h003 =  3.3285389740e-6; 
h004 = -2.8273403593e-7; 
h005 =  2.1062306160e-8; 
h006 = -2.1078768810e-9; 
h007 =  2.8019291329e-10; 
h011 = -1.5649734675e-5; 
h012 =  9.2528827145e-6; 
h013 = -3.9121289103e-7; 
h014 = -9.1317516383e-8; 
h015 =  6.2908199804e-8; 
h021 =  2.7762106484e-5; 
h022 = -5.8583034265e-6; 
h023 =  7.1016762467e-7; 
h024 =  7.1739762898e-8; 
h031 = -1.6521159259e-5; 
h032 =  3.9639828087e-6; 
h033 = -1.5377513346e-7; 
h042 = -1.7051093741e-6; 
h043 = -2.1117638838e-8; 
h041 =  6.9111322702e-6; 
h051 = -8.0539615540e-7; 
h052 =  2.5368383407e-7; 
h061 =  2.0543094268e-7; 
h101 = -3.1038981976e-4; 
h102 =  1.21312343735e-5; 
h103 = -1.9494810995e-7; 
h104 =  9.0775471288e-8; 
h105 = -2.2294250846e-8; 
h111 =  3.5009599764e-5; 
h112 = -4.7838544078e-6; 
h113 = -1.8566384852e-6; 
h114 = -6.8239240593e-8; 
h121 = -3.7435842344e-5; 
h122 = -1.18391541805e-7; 
h123 =  1.3045795693e-7; 
h131 =  2.4141479483e-5; 
h132 = -1.72793868275e-6; 
h133 =  2.5872962697e-9; 
h141 = -8.7595873154e-6; 
h142 =  6.4783588915e-7; 
h151 = -3.3052758900e-7;
h201 =  6.6928067038e-4; 
h202 = -1.7396230487e-5; 
h203 = -1.6040750532e-6; 
h204 =  4.1865759450e-9; 
h211 = -4.3592678561e-5; 
h212 =  5.5504173825e-6; 
h213 =  1.8206916278e-6; 
h221 =  3.5907822760e-5; 
h222 =  1.46416731475e-6; 
h223 = -2.1910368022e-7; 
h231 = -1.4353633048e-5; 
h232 =  1.5827653039e-7; 
h241 =  4.3703680598e-6;
h301 = -8.5047933937e-4; 
h302 =  1.87353886525e-5; 
h303 =  1.6421035666e-6; 
h311 =  3.4532461828e-5; 
h312 = -4.9223558922e-6; 
h313 = -4.5147285423e-7; 
h321 = -1.8698584187e-5; 
h322 = -2.4413069600e-7; 
h331 =  2.2863324556e-6;
h401 =  5.8086069943e-4; 
h402 = -8.6611093060e-6; 
h403 = -5.9373249090e-7; 
h411 = -1.1959409788e-5; 
h421 =  3.8595339244e-6; 
h412 =  1.2954612630e-6;
h501 = -2.1092370507e-4; 
h502 =  1.54637136265e-6; 
h511 =  1.3864594581e-6; 
h601 =  3.1932457305e-5; 

dynamic_enthalpy_part =  z.*(h001 + xs.*(h101 + xs.*(h201 + xs.*(h301 + xs.*(h401 ...
    + xs.*(h501 + h601.*xs))))) + ys.*(h011 + xs.*(h111 + xs.*(h211 + xs.*(h311 ...
    + xs.*(h411 + h511.*xs)))) + ys.*(h021 + xs.*(h121 + xs.*(h221 + xs.*(h321 ...
    + h421.*xs))) + ys.*(h031 + xs.*(h131 + xs.*(h231 + h331.*xs)) + ys.*(h041 ...
    + xs.*(h141 + h241.*xs) + ys.*(h051 + h151.*xs + h061.*ys))))) + z.*(h002 ...
    + xs.*(h102 + xs.*(h202 + xs.*(h302 + xs.*(h402 + h502.*xs)))) + ys.*(h012 ...
    + xs.*(h112 + xs.*(h212 + xs.*(h312 + h412.*xs))) + ys.*(h022 + xs.*(h122 ...
    + xs.*(h222 + h322.*xs)) + ys.*(h032 + xs.*(h132 + h232.*xs) + ys.*(h042 ...
    + h142.*xs + h052.*ys)))) + z.*(h003 + xs.*(h103 + xs.*(h203 + xs.*(h303 ...
    + h403.*xs))) + ys.*(h013 + xs.*(h113 + xs.*(h213 + h313.*xs)) + ys.*(h023 ...
    + xs.*(h123 + h223.*xs) + ys.*(h033 + h133.*xs + h043.*ys))) + z.*(h004 ...
    + xs.*(h104 + h204.*xs) + ys.*(h014 + h114.*xs + h024.*ys) + z.*(h005 ...
    + h105.*xs + h015.*ys + z.*(h006 + h007.*z))))));

dynamic_enthalpy = dynamic_enthalpy_part.*1e8;         %Note. 1e8 = db2Pa*1e4;

if transposed
    dynamic_enthalpy = dynamic_enthalpy.';
end

end
