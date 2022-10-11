function ntp_pt_vs_CT_ratio = gsw_ntp_pt_vs_CT_ratio(SA,CT,p)

% gsw_ntp_pt_vs_CT_ratio                    ratio of gradients of potential
%                             temperature and Conservative Temperature in a
%                            neutral tangent plane (in a locally-referenced
%                              potential density surface)(75-term equation)
% =========================================================================
%
% USAGE:
%  ntp_pt_vs_CT_ratio = gsw_ntp_pt_vs_CT_ratio(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the ratio of the two-dimensional gradient of potential 
%  temperature versus that of Conservative Temperature, CT, along the  
%  neutral tangent plane.  The potential temperature is the regular one  
%  which has a reference sea pressure of 0 dbar.  Part of the calculation  
%  uses the computationally-efficient 75-term expression for specific 
%  volume in terms of SA, CT and p (Roquet et al., 2015).
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
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  ntp_pt_vs_CT_ratio  =  The ratio of the spatial gradient of 
%                         potential temperature versus that of 
%                         Conservative Temperature in the 
%                         neutral tangent plane (ntp).         [ unitless ]
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
%    See Eqn. (A.14.5) of this TEOS-10 Manual.   
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_ntp_pt_vs_CT_ratio:  Requires three inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_ntp_pt_vs_CT_ratio: SA and CT must have same dimensions')
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
    error('gsw_ntp_pt_vs_CT_ratio: Inputs array dimensions arguments do not agree')
end 

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

[dummy, alpha, beta] = gsw_specvol_alpha_beta(SA,CT,p);

%--------------------------------------------------------------------------
% This function calculates the ntp_pt_vs_CT_ratio using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and 
% p.  If one wanted to compute this with the full TEOS-10 Gibbs function 
% expression for specific volume, the following lines of code will enable 
% this.
%
%    t = gsw_t_from_CT(SA,CT,p);
%    beta = gsw_beta_const_CT_t_exact(SA,t,p);
%    alpha = gsw_alpha_wrt_CT_t_exact(SA,t,p);
%
%--------- This is the end of the alternative code-------------------------

[pt_SA, pt_CT] = gsw_pt_first_derivatives(SA,CT);

ntp_pt_vs_CT_ratio = pt_CT + pt_SA.*(alpha./beta);

if transposed
    ntp_pt_vs_CT_ratio = ntp_pt_vs_CT_ratio.';
end

end
