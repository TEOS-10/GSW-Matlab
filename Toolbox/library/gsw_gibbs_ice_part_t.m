function gibbs_ice_part_t = gsw_gibbs_ice_part_t(t,p)

% gsw_gibbs_ice_part_t        part of the derivative of Gibbs energy of ice 
% =========================================================================
%
% USAGE:
%  gibbs_ice_part_t = gsw_gibbs_ice_part_t(nt,np,t,p)
%
% DESCRIPTION:
%  part of the the first temperature derivative of Gibbs energy of ice
%  that is the outout is gibbs_ice(1,0,t,p) + S0
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%
% OUTPUT:
%  gibbs_ice_part_t = part of temperature derivative       [ J kg^-1 K^-1 ]
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

rec_Pt = 1.634903221903779e-3;   % 1./Pt, where Pt = 611.657;  Experimental
% triple-point pressure in Pa.

T = t + 273.15; % The input temperature t is in-situ temperature in
                % units of degrees Celcius.  T is the in-situ Absolute
rec_Tt = 3.660858105139845e-3;   % 1/Tt = 3.660858105139845e-3; 
                % Temperature of the ice in degrees kelvin (K).
tau = T.*rec_Tt;

db2Pa = 1e4;
dzi = db2Pa.*p.*rec_Pt;

t1 = (3.68017112855051e-2 + 5.10878114959572e-2i);
t2 = (3.37315741065416e-1 + 3.35449415919309e-1i);

r1 = (4.47050716285388e1 + 6.56876847463481e1i);
r20	= (-7.25974574329220e1 - 7.81008427112870e1i);
r21	= (-5.57107698030123e-5 + 4.64578634580806e-5i);
r22	= (	2.34801409215913e-11 - 2.85651142904972e-11i);

tau_t1 = tau./t1;
tau_t2 = tau./t2;

r2 = r20 + dzi.*(r21 + r22.*dzi);

g = r1.*(log((1 + tau_t1)./(1 - tau_t1)) - 2.*tau_t1) ...
    + r2.*(log((1 + tau_t2)./(1 - tau_t2)) - 2.*tau_t2);

gibbs_ice_part_t = real(g);

end