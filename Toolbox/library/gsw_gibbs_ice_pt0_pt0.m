function gibbs_ice_pt0_pt0 = gsw_gibbs_ice_pt0_pt0(pt0)

% gsw_gibbs_ice_pt0_pt0       part of the derivative of Gibbs energy of ice 
% =========================================================================
%
% USAGE:
%  gibbs_ice_pt0_pt0 = gsw_gibbs_ice_pt0_pt0(pt0)
%
% DESCRIPTION:
%  The second temperature derivative of Gibbs energy of ice at the 
%  potential temperature with reference sea pressure of zero dbar.  That is
%  the output is gibbs_ice(2,0,pt0,0). 
%
% INPUT:
%  pt0  =  potential temperature with reference sea pressure of zero dbar
%                                                                 [ deg C ]
%
% OUTPUT:
%  gibbs_ice_pt0_pt0 = temperature second derivative at pt0
%                                                          [ J kg^-1 K^-2 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IAPWS, 2009: Revised Release on the Equation of State 2006 for H2O Ice
%   Ih. The International Association for the Properties of Water and
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

T = pt0 + 273.15;  % The input temperature is potential temperature with 
            % reference sea pressure of zero dbar units of degrees Celcius. 
% T is the Absolute potential temperature of the ice in degrees kelvin (K).
rec_Tt = 3.660858105139845e-3;   % 1/Tt = 3.660858105139845e-3; 
tau = T.*rec_Tt;

t1 = (3.68017112855051e-2 + 5.10878114959572e-2i);
t2 = (3.37315741065416e-1 + 3.35449415919309e-1i);

r1 = (4.47050716285388e1 + 6.56876847463481e1i);
r20	= (-7.25974574329220e1 - 7.81008427112870e1i);

g = r1.*(1./(t1 - tau) + 1./(t1 + tau) - 2./t1) ...
    + r20.*(1./(t2 - tau) + 1./(t2 + tau) - 2./t2);

gibbs_ice_pt0_pt0 = rec_Tt.*real(g);

end