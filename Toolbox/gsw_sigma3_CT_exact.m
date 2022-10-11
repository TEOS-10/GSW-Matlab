function sigma3_CT_exact = gsw_sigma3_CT_exact(SA,CT)

% gsw_sigma3_CT_exact                        potential density anomaly with
%                                       reference sea pressure of 3000 dbar 
%==========================================================================
% 
% USAGE:  
%  sigma3_CT_exact = gsw_sigma3_CT_exact(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 3000 
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely gsw_sigma3(SA,CT,p), 
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma3_CT_exact  =  potential density anomaly with            [ kg/m^3 ]
%                      respect to a reference pressure of 3000 dbar,   
%                      that is, this potential density - 1000 kg m^-3.
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
%    See Eqn. (A.30.1) of this TEOS-10 Manual. 
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
   error('gsw_sigma3_CT_exact:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma3_CT_exact: SA and CT must have same dimensions')
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

pr3000 = 3000*ones(size(SA));
t = gsw_t_from_CT(SA,CT,pr3000);
sigma3_CT_exact = gsw_rho_t_exact(SA,t,pr3000) - 1000;

if transposed
    sigma3_CT_exact = sigma3_CT_exact.';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% potential density with 1000 kg/m^3 subtracted from it. 

end
