function cabbeling = gsw_cabbeling(SA,CT,p)

% gsw_cabbeling                                       cabbeling coefficient 
%                                                        (48-term equation)
%==========================================================================
%
% USAGE:  
%  cabbeling = gsw_cabbeling(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the cabbeling coefficient of seawater with respect to  
%  Conservative Temperature.  This function uses the computationally-
%  efficient 48-term expression for density in terms of SA, CT and p
%  (McDougall et al., 2011)
%   
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2011).  The GSW library function 
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
%  cabbeling  =  cabbeling coefficient with respect to            [ 1/K^2 ]
%                Conservative Temperature.                    
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]   
%
% VERSION NUMBER: 3.01 (23rd March, 2011)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.9.2) and (P.4) of this TEOS-10 manual.
%
%  McDougall T.J., P.M. Barker, R. Feistel and D.R. Jackett, 2011:  A 
%   computationally efficient 48-term expression for the density of 
%   seawater in terms of Conservative Temperature, and related properties
%   of seawater.  To be submitted to Ocean Science Discussions. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_cabbeling:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_cabbeling: SA and CT must have same dimensions')
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
    error('gsw_cabbeling: Inputs array dimensions arguments do not agree')
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

dCT = 1e-3;          % increment in Conservative Temperature is 1e-3 deg C.
CT_l = CT - dCT;
CT_u = CT + dCT;

[dummy,alpha,beta] = gsw_rho_alpha_beta(SA,CT,p);
alpha_u = gsw_alpha(SA,CT_u,p);
alpha_l = gsw_alpha(SA,CT_l,p);

alpha_on_beta = alpha./beta;
alpha_CT = (alpha_u - alpha_l)./(CT_u-CT_l);

dSA = 1e-3;                % increment in Absolute Salinity is 1e-3 g kg^-1
inds_l = find(SA>=dSA);
SA_l = nan(size(SA));
if ~isempty(inds_l)
    SA_l(inds_l) = SA(inds_l) - dSA;
end
inds_l = find(SA<dSA);
if ~isempty(inds_l)
    SA_l(inds_l) = 0;
end
SA_u = SA + dSA;  

[dummy,alpha_u,beta_u] = gsw_rho_alpha_beta(SA_u,CT,p);
[dummy,alpha_l,beta_l] = gsw_rho_alpha_beta(SA_l,CT,p);

alpha_SA = (alpha_u - alpha_l)./(SA_u-SA_l);
beta_SA = (beta_u - beta_l)./(SA_u-SA_l);
cabbeling = alpha_CT + alpha_on_beta.*(2.*alpha_SA - alpha_on_beta.*beta_SA);

%--------------------------------------------------------------------------
% This function calculates cabbeling using the computationally-efficient
% 48-term expression for density in terms of SA, CT and p.  If one wanted 
% to compute cabbeling with the full TEOS-10 Gibbs function expression 
% for density, the following lines of code will enable this.
%
%   pr0 = zeros(size(p)); 
%   pt = gsw_pt_from_CT(SA,CT);    
%   dpt = 1e-3;          % increment in potential temperature is 1e-3 deg C
%   pt_l = pt - dpt; 
%   pt_u = pt + dpt; 
%   CT_l = gsw_CT_from_pt(SA,pt_l);
%   CT_u = gsw_CT_from_pt(SA,pt_u);  
%   t_l = gsw_pt_from_t(SA,pt_l,pr0,p); 
%   t_u = gsw_pt_from_t(SA,pt_u,pr0,p);
%   t = 0.5*(t_l + t_u);
%   alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p); 
%   beta = gsw_beta_const_CT_t_exact(SA,t,p); 
%   alpha_on_beta = alpha./beta;
%   alpha_CT = (gsw_alpha_wrt_CT_t_exact(SA,t_u,p)-gsw_alpha_wrt_CT_t_exact(SA,t_l,p))./(CT_u-CT_l);
%   dSA = 1e-3;                %increment in Absolute Salinity is 1e-3 g/kg
%   SA_l = nan(size(SA));
%   inds_l = find(SA>=dSA);
%   if ~isempty(inds_l)   
%     SA_l(inds_l) = SA(inds_l)-dSA;
%   end
%   inds_l = find(SA<dSA);
%   if ~isempty(inds_l)   
%     SA_l(inds_l) = 0; 
%   end
%   SA_u = SA + dSA;  
%   alpha_SA  = (gsw_alpha_wrt_CT_t_exact(SA_u,t,p)-gsw_alpha_wrt_CT_t_exact(SA_l,t,p))./(SA_u-SA_l);
%   beta_SA   = (gsw_beta_const_CT_t_exact(SA_u,t,p)-gsw_beta_const_CT_t_exact(SA_l,t,p))./(SA_u-SA_l);
%   cabbeling = alpha_CT + alpha_on_beta.*(2.*alpha_SA - alpha_on_beta.*beta_SA);
%
%---------------This is the end of the alternative code--------------------

if transposed
    cabbeling = cabbeling.';
end

end
