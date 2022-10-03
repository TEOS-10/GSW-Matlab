function [SA,CT] = gsw_SA_CT_from_sigma2_spiciness2(sigma2,spiciness2)

% gsw_SA_CT_from_sigma2_spiciness2     SA and CT from sigma2 and spiciness2
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  [SA,CT] = gsw_SA_CT_from_sigma2_spiciness2(sigma2,spiciness2)
% 
% DESCRIPTION:
%  Calculates the Absolute Salinity and the Conservative Temperature of a 
%  seawater sample at given values of potential density referenced to 
%  2000 dbar and spiciness2.  This function uses the computationally
%  efficient 75-term expression for specific volume in terms of SA, CT and
%  p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  sigma2 =  density anomaly of a seawater sample (e.g 26 kg/m^3)
%            referenced to 2000 dbar                             [ kg/m^3 ]
%   Note. This input has had 1000 kg/m^3 subtracted from it.  That is, 
%   it is 'density anomaly', not 'density'.
%
%  spiciness2 = spiciness appropriate to 2000 dbar, as given by the paper
%               of McDougall and Krzysik (2015), and the GSW algorithm 
%               gsw_spiciness2(SA,CT)                            [ kg/m^3 ]
%
%  sigma2 and spiciness2 must have the same dimensions.  
%
% OUTPUT:
%  SA  =  Absolute Salinity.                                       [ g/kg ]
%    Note, SA is expressed on the Reference-Composition 
%    Salinity Scale of Millero et al. (2008). 
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.yy (10th September 2019)
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
   error('gsw_SA_CT_from_sigma2_spiciness2:  Requires two inputs')
end 

[md,nd] = size(sigma2);
[msp,nsp] = size(spiciness2);

if (msp ~= md || nsp ~= nd)
    error('gsw_SA_CT_from_sigma2_spiciness2: sigma2 and spiciness2 must have same dimensions')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
rho = 1000.0 + sigma2; 

% initial estimate of SA from the polynomial of SA(sigma2,spiciness2)
SA = gsw_SA_poly_spiceness2(sigma2,spiciness2); 

[CT,dummy] = gsw_CT_from_rho(rho,SA,2000); % The second solution is ignored.

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure of 
% McDougall and Wotherspoon (2014)
%--------------------------------------------------------------------------
for Number_of_iterations = 1:8
    delta_spiciness = spiciness2 - gsw_spiciness2(SA,CT);
    derivative = gsw_deriv_SA_poly_spiciness2(sigma2,spiciness2);    
    SA_old = SA;
    SA = SA_old + delta_spiciness.*derivative;    
    [CT, dummy] = gsw_CT_from_rho(rho,SA,2000); % The second solution is ignored.
end

CT(SA < 0 | SA > 42) = NaN; 
SA(SA < 0 | SA > 42) = NaN;
CT(CT < -5 | CT > 40) = NaN; 
SA(CT < -5 |CT > 40) = NaN;

%--------------------------------------------------------------------------
% Note that this algorithm returns only one set of values of [SA,CT].
% At low salinities where the TMD is larger (warmer) than the freezing
% temperature there can be two solutions for the same input values 
% (sigma2,spiciness2).  
%--------------------------------------------------------------------------

end

function SA = gsw_SA_poly_spiceness2(sigma2,spiciness2)

 SAploy =[10.490055329718849
   0.633818964220827
  -0.000050331450079
   0.715503779079206
  -0.000859057970369
   0.001068697703639
  -0.000081259586352
   0.000088743890560
   0.000019621369324
  -0.000036342610189];

SA = SAploy(1) + sigma2.*(SAploy(2) + SAploy(6).*spiciness2 + sigma2.*(SAploy(3) ...
    + SAploy(7).*spiciness2 + SAploy(9)*sigma2)) + spiciness2.*(SAploy(4) ...
    + spiciness2.*(SAploy(5) + SAploy(8)*sigma2 + SAploy(10)*spiciness2));

end

function SA_poly_spiciness2_derivative = gsw_deriv_SA_poly_spiciness2(sigma2,spiciness2)

 SAploy =[10.490055329718849
   0.633818964220827
  -0.000050331450079
   0.715503779079206
  -0.000859057970369
   0.001068697703639
  -0.000081259586352
   0.000088743890560
   0.000019621369324
  -0.000036342610189];

SA_poly_spiciness2_derivative = SAploy(4) + sigma2.*(SAploy(6) + SAploy(7)*sigma2) + ...
    spiciness2.*(2.*(SAploy(5) + SAploy(8)*sigma2) + 3.*SAploy(10)*spiciness2);

end

