function pt = gsw_pt_from_CT_plus(SA,CT_plus,r)

% gsw_pt_from_CT_plus                   potential temperature from CT_plus
%==========================================================================
%
% USAGE:  
%  pt = gsw_pt_from_CT_plus(SA,CT_plus,r)
%
% DESCRIPTION:
%  Calculates potential temperature (with a reference sea pressure of
%  zero dbar) from CT_plus.  CT_plus is the version of Conservative 
%  Temperature which has an amplified difference from potential
%  temperature (pt).  That is, CT_plus obeys 
%
%        (pt - CT_plus)  =  r*(pt - CT)
%
%  This function uses NNN iterations through a modified Newton-Raphson 
%  (N-R) iterative solution proceedure, starting from a first guess as 
%  described in the Notes (paper?)  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%   r  =  amplification factor                                 [ unitless ]
%
%  SA & CT need to have the same dimensions and r should be a solitary
%  number. 
%
% OUTPUT:
%  pt  =  potential temperature referenced to a sea pressure 
%         of zero dbar (ITS-90)                                   [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker. 
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.02 (20th December, 2012)
%  It is not envisaged that this function become part of GSW.
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See sections 3.1 and 3.3 of this TEOS-10 Manual.
%
%  McDougall T.J., P.M. Barker and R. Feistel, 2013:  A computationally
%   efficient 48-term expression for the density of seawater in terms of
%   Conservative Temperature, and related properties of seawater.  
%   To be submitted to J. Atm. Ocean. Technol., xx, yyy-zzz.
%
%  McDougall T.J. and S.J. Wotherspoon, 2012: A simple modification of 
%   Newton’s method to achieve convergence of order "1 + sqrt(2)".
%   Submitted to Applied Mathematics and Computation.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_pt_from_CT_plus: Requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT_plus);

if (mt ~= ms | nt ~= ns)
    error('gsw_pt_from_CT: SA and CT must have same dimensions')
end

if ~isscalar(r)
    error('gsw_pt_from_CT: r must be a scalar')
end

if ms == 1
    SA = SA.';
    CT_plus = CT_plus.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

cp0 = 3991.86795711963;           % from Eqn. (3.3.3) of IOC et al. (2010).
SSO = 35.16504;      % this is the Reference Salinity of the Standard Ocean
a = 0.05*(1 - SA/SSO);
                   
pt = (1 - a).*CT_plus./(1 + (r - 1).*a);       % the starting estimate of pt 

df_dpt = 1 ;                              % the starting estimate of df_dpt

% start the 2 iterations through the modified Newton-Rapshon iterative 
% method (McDougall and Wotherspoon, 2012). 

for Number_of_iterations = 1:2
    pt_old = pt;
    f = pt - CT_plus -r*(pt - gsw_CT_from_pt(SA,pt));
    pt = pt_old - f./df_dpt;           % 1/2-way through a modified N-R loop
    ptm = 0.5*(pt + pt_old);
    
% This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative
% of the Gibbs function with respect to temperature at zero sea pressure.
    
    dCT_pt = -(ptm + 273.15).*gsw_gibbs_pt0_pt0(SA,ptm)./cp0;
    df_dpt = 1 - r*(1 - dCT_pt);
    pt = pt_old - f./df_dpt;           % end of a full modified N-R iteration
end

% After 2 iterations this calulation is at machine precission for this
% calculation which is 2 x 10 ^-13 (for r = 10).

if transposed
    pt = pt.';
end


end
