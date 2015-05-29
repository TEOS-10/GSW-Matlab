function thermobaric_CT_exact = gsw_thermobaric_CT_exact(SA,CT,p)

% gsw_thermobaric_CT_exact                          thermobaric coefficient 
%==========================================================================
%
% USAGE:  
%  thermobaric_CT_exact = gsw_thermobaric_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates the thermobaric
%  coefficient with the full TEOS-10 Gibbs function expression for density.
%  This function uses finite differences to calculate the temperature and
%  pressure derivatives.
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely gsw_thermobaric(SA,CT,p) 
%  which uses the computationally efficient 48-term expression for density 
%  in terms of SA, CT and p (IOC et al., 2010).   
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
%  thermobaric_CT_exact  =  thermobaric coefficient with       [ 1/(K Pa) ] 
%                           respect to Conservative Temperature.           
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (5th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_thermobaric_CT_exact: Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_thermobaric_CT_exact: SA and CT must have same dimensions')
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
    error('gsw_thermobaric_CT_exact: Inputs array dimensions arguments do not agree')
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

beta = gsw_beta_CT_exact(SA,CT,p);

dp = 0.1;                        % pressure increment is 0.1 dbar (1000 Pa)
p_u = p - dp;
p_l = p + dp;

alpha_on_beta_u = gsw_alpha_on_beta_CT_exact(SA,CT,p_u);
alpha_on_beta_l = gsw_alpha_on_beta_CT_exact(SA,CT,p_l);

thermobaric_CT_exact = beta.*(alpha_on_beta_u - alpha_on_beta_l)./(p_u - p_l);
Pa2db = 1e-4;  % To have units of 1/(K Pa)
thermobaric_CT_exact = thermobaric_CT_exact.*Pa2db;  


if transposed
    thermobaric_CT_exact = thermobaric_CT_exact.';
end

end
