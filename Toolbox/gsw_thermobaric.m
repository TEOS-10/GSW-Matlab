function thermobaric = gsw_thermobaric(SA,CT,p)

% gsw_thermobaric                thermobaric coefficient (48-term equation)
%==========================================================================
%
% USAGE:  
%  thermobaric = gsw_thermobaric(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates rho from the 
%  computationally-efficient 48-term expression for density in terms of
%  SA, CT and p (McDougall et al., 2013).
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2013).  The GSW library function 
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
%  thermobaric  =  thermobaric coefficient with                [ 1/(K Pa) ] 
%                  respect to Conservative Temperature.           
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.02 (16th November, 2012)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
%
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2013:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to J. Atm. Ocean. Technol., xx, yyy-zzz.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_thermobaric: Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_thermobaric: SA and CT must have same dimensions')
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
    error('gsw_thermobaric: Inputs array dimensions arguments do not agree')
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

db2Pa = 1e4;
dp = 1e-1;                        % pressure increment is 1e-1 dbar (10 Pa)

p_u = zeros(size(p));
p_u(p >= dp) = p(p >= dp) - dp;

p_l = dp*ones(size(p));
p_l(p >= dp) = p(p >= dp) + dp;

[dummy,alpha,beta] = gsw_rho_alpha_beta(SA,CT,p);
[dummy,alpha_u,beta_u] = gsw_rho_alpha_beta(SA,CT,p_u);
[dummy,alpha_l,beta_l] = gsw_rho_alpha_beta(SA,CT,p_l);

alpha_p = (alpha_u - alpha_l)./(p_u - p_l);
beta_p  = (beta_u - beta_l)./(p_u - p_l);

thermobaric = alpha_p - (alpha./beta).*beta_p;
thermobaric = thermobaric./db2Pa;         % To have units of 1/(K Pa)

%--------------------------------------------------------------------------
% This function calculates thermobaric_CT using the computationally
% efficient 48-term expression for density in terms of SA, CT and p.  If 
% one wanted to compute thermobaric_CT with the full TEOS-10 Gibbs function
% expression for density, the following lines of code will do this.
%
%  pr0 = zeros(size(p)); 
%  pt = gsw_pt_from_CT(SA,CT);
%  t_l = gsw_pt_from_t(SA,pt,pr0,p_l);   
%  t_u = gsw_pt_from_t(SA,pt,pr0,p_u);
%  t = 0.5*(t_l + t_u);
%  alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%  beta = gsw_beta_const_CT_t_exact(SA,t,p);
%  alpha_p = (gsw_alpha_wrt_CT_t_exact(SA,t_u,p_u) - gsw_alpha_wrt_CT_t_exact(SA,t_l,p_l))./(p_u - p_l);
%  beta_p = (gsw_beta_const_CT_t_exact(SA,t_u,p_u) - gsw_beta_const_CT_t_exact(SA,t_l,p_l))./(p_u - p_l);
%  thermobaric = alpha_p - (alpha./beta).*beta_p;
%  thermobaric = thermobaric./db2Pa;      % To have units of 1/(K Pa)
%
%----------------This is the end of the alternative code-------------------

if transposed
    thermobaric = thermobaric.';
end

end
