function [IPV_vs_fNsquared_ratio, p_mid] = gsw_IPV_vs_fNsquared_ratio(SA,CT,p,p_ref)

% gsw_IPV_vs_fNsquared_ratio              ratio of the vertical gradient of
%                       potential density (with reference pressure, p_ref),  
%                            to the vertical gradient of locally-referenced 
%                                      potential density (75-term equation)
%==========================================================================
% 
% USAGE:  
%  [IPV_vs_fNsquared_ratio, p_mid] = gsw_IPV_vs_fNsquared_ratio(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates the ratio of the vertical gradient of potential density to 
%  the vertical gradient of locally-referenced potential density.  This 
%  ratio is also the ratio of the planetary Isopycnal Potential Vorticity
%  (IPV) to f times N^2, hence the name for this variable,
%  IPV_vs_fNsquared_ratio (see Eqn. (3.20.17) of IOC et al. (2010)). 
%  The reference sea pressure, p_ref, of the potential density surface must
%  have a constant value.
%
%  IPV_vs_fNsquared_ratio is evaluated at the mid pressure between the 
%  individual data points in the vertical.  This function uses the 
%  computationally-efficient 75-term expression for specific volume in 
%  terms of SA, CT and p (Roquet et al., 2015). 
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:  
%  SA    = Absolute Salinity                                       [ g/kg ]
%  CT    = Conservative Temperature (ITS-90)                      [ deg C ]
%  p     = sea pressure                                            [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref = reference sea pressure of the potential density surface
%         ( i.e. absolute reference pressure - 10.1325 dbar )      [ dbar ]
%
%  SA & CT need to have the same dimensions.
%  p & p_ref may have dimensions 1x1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  IPV_vs_fNsquared_ratio = The ratio of the vertical gradient of 
%          potential density referenced to p_ref, to the vertical 
%          gradient of locally-referenced potential density.  It is 
%          output on the same vertical (M-1)xN grid as p_mid. 
%          IPV_vs_fNsquared_ratio is dimensionless.            [ unitless ]
%  p_mid = mid pressure between the individual points of the p grid.
%          That is, p_mid is on a (M-1)xN grid.                    [ dbar ]
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
%    See Eqn. (3.20.5) of this TEOS-10 Manual. 
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
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 | nargin == 4)
   error('gsw_IPV_vs_fNsquared_ratio:  Requires three or four inputs')
end %if
if ~(nargout == 2)
   error('gsw_IPV_vs_fNsquared_ratio:  Requires two outputs')
end %if

if nargin == 3
%    Assume reference pressure, p_ref, is 0 dbar.
  p_ref = 0;
end %if

if ~isscalar(unique(p_ref))
    error('gsw_IPV_vs_fNsquared_ratio: The reference pressures differ, they should be unique')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_IPV_vs_fNsquared_ratio: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('gsw_IPV_vs_fNsquared_ratio: There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p scalar - must be two bottles
    error('gsw_IPV_vs_fNsquared_ratio: There must be at least 2 values')
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
    error('gsw_IPV_vs_fNsquared_ratio: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    [mp,np] = size(p);
    transposed = 1;
else
    transposed = 0;
end

p_ref = unique(p_ref)*ones(mp-1,np);               %resize the reference pressure 

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

Ishallow = 1:(mp-1);
Ideep = 2:mp;
p_mid = 0.5*(p(Ishallow,:) + p(Ideep,:));
SA_mid = 0.5*(SA(Ishallow,:) + SA(Ideep,:));
CT_mid = 0.5*(CT(Ishallow,:) + CT(Ideep,:));

dSA = SA(Ishallow,:) - SA(Ideep,:);
dCT = CT(Ishallow,:) - CT(Ideep,:);

[dummy,alpha,beta] = gsw_rho_alpha_beta(SA_mid,CT_mid,p_mid);
[dummy,alpha_pref,beta_pref] = gsw_rho_alpha_beta(SA_mid,CT_mid,p_ref);

%--------------------------------------------------------------------------
% This function calculates IPV_vs_fNsquared_ratio using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and 
% p.  If one wanted to compute this with the full TEOS-10 Gibbs function
% expression for specific volume, the following lines of code will enable
% this.
%
%    t_mid = gsw_pt_from_CT(SA_mid,CT_mid,p_mid);
%    beta = gsw_beta_const_CT_t_exact(SA_mid,t_mid,p_mid);
%    alpha  = gsw_alpha_wrt_CT_t_exact(SA_mid,t_mid,p_mid);
%    beta_pref  = gsw_beta_const_CT_t_exact(SA_mid,t_mid,p_ref);
%    alpha_pref = gsw_alpha_wrt_CT_t_exact(SA_mid,t_mid,p_ref);
%
%-----------This is the end of the alternative code------------------------ 

numerator = dCT.*alpha_pref - dSA.*beta_pref;
denominator = dCT.*alpha - dSA.*beta;

IPV_vs_fNsquared_ratio = nan(size(SA_mid));
IPV_vs_fNsquared_ratio(denominator ~= 0) = numerator(denominator ~= 0)./denominator(denominator ~= 0);

if transposed
    IPV_vs_fNsquared_ratio = IPV_vs_fNsquared_ratio.';
    p_mid = p_mid.';
end

end