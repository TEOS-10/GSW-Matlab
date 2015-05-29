function sigma1_CT = gsw_sigma1_CT(SA,CT)

% gsw_sigma1_CT                              potential density anomaly with
%                                      reference sea pressure of 1000 dbar.
%==========================================================================
% 
% USAGE:  
%  sigma1_CT = gsw_sigma1_CT(SA,CT)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference pressure of 1000 
%  dbar, this being this particular potential density minus 1000 kg/m^3.
%  This function has inputs of Absolute Salinity and Conservative
%  Temperature.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  sigma1_CT  =  potential density anomaly with                  [ kg/m^3 ]
%                respect to a reference pressure of 1000 dbar,   
%                that is, this potential density - 1000 kg/m^3.
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
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
   error('gsw_sigma1_CT:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma1_CT: SA and CT must have same dimensions')
end

if ms == 1
    SA = SA';
    CT = CT';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pr0 = zeros(size(SA));
pr1000 = 1000*ones(size(SA));
pt0 = gsw_pt_from_CT(SA,CT);
t1000 = gsw_pt_from_t(SA,pt0,pr0,pr1000);

n0 = 0;
n1 = 1;

sigma1_CT = ones(size(SA))./gsw_gibbs(n0,n0,n1,SA,t1000,pr1000) ...
               - 1000;

if transposed
    sigma1_CT = sigma1_CT';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% potential density with 1000 kg/m^3 subtracted from it. 

end
