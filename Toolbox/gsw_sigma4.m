function sigma4 = gsw_sigma4(SA,CT)

% gsw_sigma4                       potential density anomaly with reference
%                              sea pressure of 4000 dbar (75-term equation)
%==========================================================================
% 
% USAGE:  
%  sigma4 = gsw_sigma4(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 4000 
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  Temperature.  This function uses the computationally-efficient 
%  expression for specific volume in terms of SA, CT and p (Roquet et al.,
%  2015).
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
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma4  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 4000 dbar,   
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual. 
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
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_sigma4:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma4: SA and CT must have same dimensions')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pr4000 = 4000*ones(size(SA));

rho4 = gsw_rho(SA,CT,pr4000);

%--------------------------------------------------------------------------
% This function calculates rho using the computationally-efficient 
% 75-term expression for specific volume in terms of SA, CT and p.  If one
% wanted to compute rho with the full TEOS-10 Gibbs function expression for 
% specific volume, the following lines of code will enable this.
%
%  rho4 = gsw_rho_CT_exact(SA,CT,pr4000);
%
%---------------This is the end of the alternative code -------------------

sigma4 = rho4 - 1000;

if transposed
    sigma4 = sigma4.';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% potential density with 1000 kg/m^3 subtracted from it. 

end
