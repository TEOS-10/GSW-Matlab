function sigma2 = gsw_sigma2(SA,CT)

% gsw_sigma2                       potential density anomaly with reference
%                              sea pressure of 2000 dbar (48-term equation)
%==========================================================================
% 
% USAGE:  
%  sigma2 = gsw_sigma2(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 2000 
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that the 48-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in IOC et al. (2010).  The GSW library function 
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
%  sigma2  =  potential density anomaly with                     [ kg/m^3 ]
%             respect to a reference pressure of 2000 dbar,   
%             that is, this potential density - 1000 kg/m^3.
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (A.30.1) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_sigma2:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma2: SA and CT must have same dimensions')
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

pr2000 = 2000*ones(size(SA));

rho2 = gsw_rho(SA,CT,pr2000);

%--------------------------------------------------------------------------
% This function calculates rho using the computationally-efficient 
% 48-term expression for density in terms of SA, CT and p.  If one wanted 
% to compute rho with the full TEOS-10 Gibbs function expression for 
% density, the following lines of code will enable this.
%
%  rho2 = gsw_rho_CT_exact(SA,CT,pr2000);
%
%---------------This is the end of the alternative code -------------------

sigma2 = rho2 - 1000;

if transposed
    sigma2 = sigma2.';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% potential density with 1000 kg/m^3 subtracted from it. 

end
