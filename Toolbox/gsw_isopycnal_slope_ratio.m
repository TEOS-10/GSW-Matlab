function isopycnal_slope_ratio = gsw_isopycnal_slope_ratio(SA,CT,p,p_ref)

% gsw_isopycnal_slope_ratio               ratio of the slopes of isopycnals
%                                      on the SA-CT diagram for p and p_ref
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
% isopycnal_slope_ratio = gsw_isopycnal_slope_ratio(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the ratio of alpha/beta at pressure, p, to that at reference
%  pressure, p_ref.  This function uses the computationally-efficient 
%  75-term expression for specific volume in terms of SA, CT and p 
%  (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of 
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
%  pr  =  reference pressure                                       [ dbar ]
%         ( i.e. absolute reference pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p and p_ref may have dimensions 1x1 or Mx1 or 1xN or MxN, where 
%  SA and CT are MxN
%
% OUTPUT:
%  isopycnal_slope_ratio  
%               =  The ratio of alpha/beta evaluated at        [ unitless ]
%                  pressure, p, to that at reference pressure, p_ref.        
%
% AUTHOR: 
%  Trevor McDougall, Paul Barker & David Jackett       [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqn. (3.17.2) of this TEOS-10 Manual.   
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
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_isopycnal_slope_ratio:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mpr,npr] = size(p_ref);

if (mt ~= ms | nt ~= ns)
    error('gsw_isopycnal_slope_ratio: SA and CT must have same dimensions')
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
    error('gsw_isopycnal_slope_ratio: Inputs array dimensions arguments do not agree')
end %if

if (mpr == 1) & (npr == 1)              % p_ref scalar - fill to size of SA
    p_ref = p_ref*ones(size(SA));
elseif (ns == npr) & (mpr == 1)         % p_ref is row vector,
    p_ref = p_ref(ones(1,ms), :);              % copy down each column.
elseif (ms == mpr) & (npr == 1)         % p_ref is column vector,
    p_ref = p_ref(:,ones(1,ns));               % copy across each row.
elseif (ns == mpr) & (npr == 1)          % p_ref is a transposed row vector,
    p_ref = p_ref.';                              % transposed then
    p_ref = p_ref(ones(1,ms), :);                % copy down each column.
elseif (ms == mpr) & (ns == npr)
    % ok
else
    error('gsw_isopycnal_slope_ratio: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    p_ref = p_ref.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

[dummy, alpha, beta] = gsw_specvol_alpha_beta(SA,CT,p);
[dummy, alpha_pref, beta_pref] = gsw_specvol_alpha_beta(SA,CT,p_ref);

%--------------------------------------------------------------------------
% This function calculates isopycnal_slope_ratio using the computationally
% efficient 75-term expression for specific volume as a function of SA, CT
% and p.  If one wanted to compute this with the full TEOS-10 Gibbs 
% function expression for specific volume, the following lines of code will 
% enable this.
%  
%     t = gsw_pt_from_CT(SA,CT,p);
%     alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%     beta = gsw_beta_const_CT_t_exact(SA,t,p);
%     tr = gsw_pt_from_t(SA,pt,p_ref0,p_ref);
%     alpha_pref = gsw_alpha_wrt_CT_t_exact(SA,tr,p_ref);
%     beta_pref = gsw_beta_const_CT_t_exact(SA,tr,p_ref);
%
%--------------This is the end of the alternative code---------------------

isopycnal_slope_ratio = NaN(size(SA));
isopycnal_slope_ratio(alpha_pref ~= 0) = (alpha(alpha_pref ~= 0).*beta_pref(alpha_pref ~= 0))./ ...
                                           (alpha_pref(alpha_pref ~= 0).*beta(alpha_pref ~= 0));

if transposed
    isopycnal_slope_ratio = isopycnal_slope_ratio.';
end

end
