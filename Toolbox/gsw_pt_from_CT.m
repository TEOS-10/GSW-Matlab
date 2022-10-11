function pt = gsw_pt_from_CT(SA,CT)

% gsw_pt_from_CT        potential temperature from Conservative Temperature
%==========================================================================
%
% USAGE:  
%  pt = gsw_pt_from_CT(SA,CT)
%
% DESCRIPTION:
%  Calculates potential temperature (with a reference sea pressure of
%  zero dbar) from Conservative Temperature.  This function uses 1.5 
%  iterations through a modified Newton-Raphson (N-R) iterative solution 
%  proceedure, starting from a rational-function-based initial condition 
%  for both pt and dCT_dpt. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt  =  potential temperature referenced to a sea pressure 
%         of zero dbar (ITS-90)                                   [ deg C ]
%
% AUTHOR: 
%  Trevor McDougall, David Jackett, Claire Roberts-Thomson and Paul Barker. 
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See sections 3.1 and 3.3 of this TEOS-10 Manual.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_pt_from_CT: Requires two inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_pt_from_CT: SA and CT must have same dimensions')
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

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

s1 = SA*0.995306702338459;   % Note that 0.995306702338459 = (35./35.16504) 

a0 = -1.446013646344788e-2;     
a1 = -3.305308995852924e-3;     
a2 =  1.062415929128982e-4;     
a3 =  9.477566673794488e-1;    
a4 =  2.166591947736613e-3;
a5 =  3.828842955039902e-3;

b0 =  1;
b1 =  6.506097115635800e-4;
b2 =  3.830289486850898e-3;
b3 =  1.247811760368034e-6;

a5CT = a5*CT; 
b3CT = b3*CT;
CT_factor = (a3 + a4*s1 + a5CT);
pt_num = a0 + s1.*(a1 + a2*s1) + CT.*CT_factor;
pt_recden = 1./(b0 + b1*s1 + CT.*(b2 + b3CT));
pt = pt_num.*pt_recden;
                          % At this point the abs max error is 1.5e-2 deg C

dpt_dCT = (CT_factor + a5CT - (b2 + b3CT + b3CT).*pt).*pt_recden;
                          
% start the 1.5 iterations through the modified Newton-Rapshon iterative 
% method (McDougall and Wotherspoon, 2014). 

CT_diff = gsw_CT_from_pt(SA,pt) - CT;
pt_old = pt;
pt = pt_old - CT_diff.*dpt_dCT; % 1/2-way through the 1st modified N-R loop
                          % At this point the abs max error is 6.6e-5 deg C

ptm = 0.5*(pt + pt_old);

% This routine calls gibbs_pt0_pt0(SA,pt0) to get the second derivative 
% of the Gibbs function with respect to temperature at zero sea pressure.  
dpt_dCT = -gsw_cp0./((ptm + gsw_T0).*gsw_gibbs_pt0_pt0(SA,ptm));
pt = pt_old - CT_diff.*dpt_dCT;  % end of 1st full modified N-R iteration
                          % At this point the abs max error is 1.e-10 deg C

CT_diff = gsw_CT_from_pt(SA,pt) - CT;
pt_old = pt;
pt = pt_old - CT_diff.*dpt_dCT; % 1.5 iterations of the modified N-R method

if transposed
    pt = pt.';
end

% The abs max error of the result is 1.42e-14 deg C

end
