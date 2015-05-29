function pt = gsw_pt_from_entropy(SA,entropy)

% gsw_pt_from_entropy                          potential temperature with a
%                                       reference sea pressure of zero dbar
%                                                  as a function of entropy 
% =========================================================================
%
% USAGE:
%  pt = gsw_pt_from_entropy(SA,entropy)
%
% DESCRIPTION:
%  Calculates potential temperature with reference pressure p_ref = 0 dbar 
%  and with entropy as an input variable. 
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  entropy  =  specific entropy                                   [ deg C ]
%
%  SA & entropy need to have the same dimensions.
%
% OUTPUT:
%  pt   =  potential temperature                                  [ deg C ]
%          with reference sea pressure (p_ref) = 0 dbar.
%  Note. The reference sea pressure of the output, pt, is zero dbar.
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix  A.10 of this TEOS-10 Manual. 
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

if ~(nargin == 2)
    error('gsw_pt_from_entropy: Requires 2 inputs, Absolute Salinity and entropy')
end %if

[ms,ns] = size(SA);
[me,ne]   = size(entropy);

if (ms ~= me | ns ~= ne )
    error('gsw_pt_from_entropy: Input arguments do not have the same dimensions')
end %if

if ms == 1
    SA = SA.';
    entropy = entropy.';
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
SSO = 35.16504;                    % from section 2.4 of IOC et al. (2010).

%  Find the initial value of pt
part1 = 1 - SA./SSO;
part2 = 1 - 0.05.*part1;
ent_SA = (cp0/273.15).*part1.*(1 - 1.01.*part1);
c = (entropy - ent_SA).*part2./cp0;
pt = 273.15*(exp(c) - 1);
dentropy_dt = cp0./((273.15 + pt).*part2); %this is the intial value of dentropy_dt

for Number_of_iterations = 1:2
    pt_old = pt;
    dentropy = gsw_entropy_from_pt(SA,pt_old) - entropy;
    pt = pt_old - dentropy./dentropy_dt ; % this is half way through the modified method (McDougall and Wotherspoon, 2012)
    ptm = 0.5*(pt + pt_old);
    dentropy_dt = -gsw_gibbs_pt0_pt0(SA,ptm);
    pt = pt_old - dentropy./dentropy_dt;
end
    
if transposed
    pt = pt.';
end

% maximum error of 2.2x10^-6 degrees C for one iteration.
% maximum error is 1.4x10^-14 degrees C for two iterations 
% (two iterations is the default, "for Number_of_iterations = 1:2"). 

end
