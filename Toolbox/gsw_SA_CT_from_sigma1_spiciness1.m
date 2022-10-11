function [SA,CT] = gsw_SA_CT_from_sigma1_spiciness1(sigma1,spiciness1)

% gsw_SA_CT_from_sigma1_spiciness1     SA and CT from sigma1 and spiciness1
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  [SA,CT] = gsw_SA_CT_from_sigma1_spiciness1(sigma1,spiciness1)
% 
% DESCRIPTION:
%  Calculates the Absolute Salinity and the Conservative Temperature of a 
%  seawater sample at given values of potential density referenced to 
%  1000 dbar and spiciness1.  This function uses the computationally
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
%  sigma1 =  density anomaly of a seawater sample (e.g 26 kg/m^3)
%            referenced to 1000 dbar                             [ kg/m^3 ]
%   Note. This input has had 1000 kg/m^3 subtracted from it.  That is, 
%   it is 'density anomaly', not 'density'.
%
%  spiciness1 = spiciness appropriate to 1000 dbar, as given by the paper
%               of McDougall and Krzysik (2015), and the GSW algorithm 
%               gsw_spiciness1(SA,CT)                            [ kg/m^3 ]
%
%  sigma1 and spiciness1 must have the same dimensions.  
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
   error('gsw_SA_CT_from_sigma1_spiciness1:  Requires two inputs')
end 

[md,nd] = size(sigma1);
[msp,nsp] = size(spiciness1);

if (msp ~= md || nsp ~= nd)
    error('gsw_SA_CT_from_sigma1_spiciness1: sigma1 and spiciness1 must have same dimensions')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
rho = 1000.0 + sigma1; 

% initial estimate of SA from the polynomial of SA(sigma1,spiciness1)
SA = gsw_SA_poly_spiceness1(sigma1,spiciness1); 

[CT,dummy] = gsw_CT_from_rho(rho,SA,1000); % The second solution is ignored.

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson iterative procedure of 
% McDougall and Wotherspoon (2014)
%--------------------------------------------------------------------------
for Number_of_iterations = 1:7
    delta_spiciness = spiciness1 - gsw_spiciness1(SA,CT);   
    derivative = gsw_deriv_SA_poly_spiciness1(sigma1,spiciness1);    
    SA_old = SA;
    SA = SA_old + delta_spiciness.*derivative;
    [CT, dummy] = gsw_CT_from_rho(rho,SA,1000); % The second solution is ignored. 
end

CT(SA < 0 | SA > 42) = NaN; 
SA(SA < 0 | SA > 42) = NaN;
CT(CT < -5 | CT > 40) = NaN; 
SA(CT < -5 |CT > 40) = NaN;

%--------------------------------------------------------------------------
% Note that this algorithm returns only one set of values of [SA,CT].
% At low salinities where the TMD is larger (warmer) than the freezing
% temperature there can be two solutions for the same input values 
% (sigma1,spiciness1).  
%--------------------------------------------------------------------------

end

function SA = gsw_SA_poly_spiceness1(sigma1,spiciness1)

 SApoly =[13.695625022104206
   0.627353321828843
   0.000074159905817
   0.711767497207624
  -0.000682109830188
   0.000726104526580
  -0.000091267411622
   0.000099817118989
   0.000022649012391
  -0.000039992686627];
 
SA = SApoly(1) + sigma1.*(SApoly(2) + SApoly(6).*spiciness1 + sigma1.*(SApoly(3) ...
    + SApoly(7).*spiciness1 + SApoly(9)*sigma1)) + spiciness1.*(SApoly(4) ...
    + spiciness1.*(SApoly(5) + SApoly(8)*sigma1 + SApoly(10)*spiciness1));

end

function SA_poly_spiciness1_derivative = gsw_deriv_SA_poly_spiciness1(sigma1,spiciness1)

 SApoly =[13.695625022104206
   0.627353321828843
   0.000074159905817
   0.711767497207624
  -0.000682109830188
   0.000726104526580
  -0.000091267411622
   0.000099817118989
   0.000022649012391
  -0.000039992686627];

 SA_poly_spiciness1_derivative = SApoly(4) + sigma1.*(SApoly(6) + SApoly(7)*sigma1) + ...
    spiciness1.*(2.*(SApoly(5) + SApoly(8)*sigma1) + 3.*SApoly(10)*spiciness1);

end

