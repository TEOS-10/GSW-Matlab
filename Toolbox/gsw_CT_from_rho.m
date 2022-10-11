function [CT,CT_multiple] = gsw_CT_from_rho(rho,SA,p)

% gsw_CT_from_rho                     Conservative Temperature from density
%                                                        (75-term equation)
% =========================================================================
%
% USAGE:
%  [CT,CT_multiple] = gsw_CT_from_rho(rho,SA,p)
%
% DESCRIPTION:
%  Calculates the Conservative Temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar), 
%  using the computationally-efficient expression for specific volume in 
%  terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  CT  =  Conservative Temperature  (ITS-90)                      [ deg C ]
%  CT_multiple  =  Conservative Temperature  (ITS-90)             [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      Conservative Temperatures for a single density.  This programme will
%      output both valid solutions.  To see this second solution the user 
%      must call the programme with two outputs (i.e. [CT,CT_multiple]), if
%      there is only one possible solution and the programme has been 
%      called with two outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
    error('gsw_CT_from_rho:  Requires three inputs')
end %if

[md,nd] = size(rho);
[ms,ns] = size(SA);
[mp,np] = size(p);

if (ms ~= md || ns ~= nd)
    error('gsw_CT_from_rho: rho and SA must have same dimensions')
end

if (mp == 1) & (np == 1)                   % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)              % p is row vector,
    p = p(ones(1,md), :);                  % copy down each column.
elseif (md == mp) & (np == 1)              % p is column vector,
    p = p(:,ones(1,nd));                   % copy across each row.
elseif (nd == mp) & (np == 1)              % p is a transposed row vector,
    p = p.';                               % transposed then
    p = p(ones(1,md), :);                  % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
   error('gsw_CT_from_rho: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA<0 | SA>42 | p <-1.5 | p>12000) = NaN;   % SA out of range, set to NaN

rho_40 = gsw_rho(SA,40*ones(size(SA)),p);
SA((rho - rho_40) < 0) = NaN;                   % rho too light, set to NaN

if any(isnan(SA(:) + p(:) + rho(:)))
    [I_bad] = find(isnan(SA + p + rho));
    SA(I_bad) = NaN;
end

CT_max_rho = gsw_CT_maxdensity(SA,p);
max_rho = gsw_rho(SA,CT_max_rho,p);
[bla,bla,rho_CT_CT,bla,bla] = gsw_rho_second_derivatives(SA,CT_max_rho,p);
deriv_lin = (rho_40 - max_rho)./(40 - CT_max_rho);
discrim = deriv_lin.*deriv_lin - 4.*0.5.*rho_CT_CT.*(max_rho - rho);
discrim(discrim < 0) = 0;
CT = CT_max_rho + 2.*(max_rho - rho)./(-deriv_lin + sqrt(discrim)); 

%    Having found this initial value of CT, begin the iterative solution
%    for the CT_upper part, the solution warmer than the TMD
for Number_of_iterations = 1:4 
 [dummy,rho_CT,dummy] = gsw_rho_first_derivatives(SA,CT,p);
 [dummy,dummy,rho_CT_CT,dummy,dummy] = gsw_rho_second_derivatives(SA,CT,p);
      b = rho_CT;
      a = 0.5.*rho_CT_CT;
      c = gsw_rho(SA,CT,p) - rho; 
      discrim = b.*b - 4.*a.*c;
      discrim(discrim < 0) = 0;
  CT_old = CT; 
  CT = CT_old + (2.*c)./(-b + sqrt(discrim));
end
  CT_upper = CT; 
   
%    Now start the CT_multiple part, the solution cooler than the TMD
CT = 2.*CT_max_rho - CT_upper;
for Number_of_iterations = 1:3                     
 [dummy,rho_CT,dummy] = gsw_rho_first_derivatives(SA,CT,p);
 [dummy,dummy,rho_CT_CT,dummy,dummy] = gsw_rho_second_derivatives(SA,CT,p); 
      b = rho_CT;
      a = 0.5.*rho_CT_CT;
      c = gsw_rho(SA,CT,p) - rho;  
      discrim = b.*b - 4.*a.*c;
      discrim(discrim < 0) = 0;
  CT_old = CT; 
  CT = CT_old + (2.*c)./(-b - sqrt(discrim)); 
                                    % Note the sign change of the sqrt term  
end
  CT_multiple = CT; 
    
% set values outside the relevant bounds to NaNs 
CT_freezing = gsw_CT_freezing(SA,p,0);  % This assumes that the seawater is 
                                        % always unsaturated with air
                                       
CT_upper(CT_upper < CT_freezing) = NaN;
CT_upper(CT_upper < CT_max_rho) = NaN;
CT_multiple(CT_multiple < CT_freezing) = NaN;
CT_multiple(CT_multiple > CT_max_rho) = NaN;

   CT = CT_upper;

if transposed
    CT = CT.';
    CT_multiple = CT_multiple.';
end

end
