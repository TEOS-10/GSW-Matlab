function [SA,CT] = gsw_SA_CT_from_sigma0_spiciness0(sigma0,spiciness0)

% gsw_SA_CT_from_sigma0_spiciness0     SA and CT from sigma0 and spiciness0
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  [SA,CT] = gsw_SA_CT_from_sigma0_spiciness0(sigma0,spiciness0)
% 
% DESCRIPTION:
%  Calculates the Absolute Salinity and the Conservative Temperature of a 
%  seawater sample at given values of potential density referenced to 
%  0 dbar and spiciness0.  This function uses the computationally-efficient 
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
%  sigma0 =  density anomaly of a seawater sample (e.g 26 kg/m^3)
%            referenced to 0 dbar                                [ kg/m^3 ]
%   Note. This input has had 1000 kg/m^3 subtracted from it.  That is, 
%   it is 'density anomaly', not 'density'.
%
%  spiciness0 = spiciness appropriate to 0 dbar, as given by the paper of
%               McDougall and Krzysik (2015), and the GSW algorithm 
%               gsw_spiciness0(SA,CT)                            [ kg/m^3 ]
%
%  sigma0 and spiciness0 must have the same dimensions.  
%
% OUTPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (2nd July 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., and O.A. Krzysik, 2015: Spiciness. Journal of Marine 
%   Research, 73, 141-152.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==2)
   error('gsw_SA_CT_from_sigma0_spiciness0:  Requires two inputs')
end 

[md,nd] = size(sigma0);
[msp,nsp] = size(spiciness0);

if (msp ~= md || nsp ~= nd)
    error('gsw_SA_CT_from_sigma0_spiciness0: sigma0 and spiciness0 must have same dimensions')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
rho = 1000.0 + sigma0; 

% initial estimate of SA from the polynomial of SA(sigma0,spiciness0)
SA = gsw_SA_poly_spiceness0(sigma0,spiciness0); 

[CT,dummy] = gsw_CT_from_rho(rho,SA,0); % The second solution is ignored.

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure of 
% McDougall and Wotherspoon (2014)
%--------------------------------------------------------------------------
for Number_of_iterations = 1:7
    delta_spiciness = spiciness0 - gsw_spiciness0(SA,CT);   
    derivative = gsw_deriv_SA_poly_spiciness0(sigma0,spiciness0);    
    SA_old = SA;
    SA = SA_old + delta_spiciness.*derivative;    
    [CT, dummy] = gsw_CT_from_rho(rho,SA,0); % The second solution is ignored.
end

CT(SA < 0 | SA > 42) = NaN; 
SA(SA < 0 | SA > 42) = NaN;
CT(CT < -5 | CT > 40) = NaN; 
SA(CT < -5 |CT > 40) = NaN;

%--------------------------------------------------------------------------
% Note that this algorithm returns only one set of values of [SA,CT].
% At low salinities where the TMD is larger (warmer) than the freezing
% temperature there can be two solutions for the same input values 
% (sigma0,spiciness0).  
%--------------------------------------------------------------------------

end

function SA = gsw_SA_poly_spiceness0(sigma0,spiciness0)

SAploy =[16.907145985921161
   0.622223077618381
   0.000223330499353
   0.706643311920640
  -0.000454964160674
   0.000312654210024
  -0.000103973707417
   0.000113069609374
   0.000026554807951
  -0.000044377349558];

SA = SAploy(1) + sigma0.*(SAploy(2) + SAploy(6).*spiciness0 + sigma0.*(SAploy(3) ...
    + SAploy(7).*spiciness0 + SAploy(9)*sigma0)) + spiciness0.*(SAploy(4) ...
    + spiciness0.*(SAploy(5) + SAploy(8)*sigma0 + SAploy(10)*spiciness0));

end

function SA_poly_spiciness0_derivative = gsw_deriv_SA_poly_spiciness0(sigma0,spiciness0)

SAploy =[16.907145985921161
   0.622223077618381
   0.000223330499353
   0.706643311920640
  -0.000454964160674
   0.000312654210024
  -0.000103973707417
   0.000113069609374
   0.000026554807951
  -0.000044377349558];
 
SA_poly_spiciness0_derivative = SAploy(4) + sigma0.*(SAploy(6) + SAploy(7)*sigma0) + ...
    spiciness0.*(2.*(SAploy(5) + SAploy(8)*sigma0) + 3.*SAploy(10)*spiciness0);

end

