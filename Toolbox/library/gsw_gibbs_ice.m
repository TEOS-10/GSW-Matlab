function gibbs_ice = gsw_gibbs_ice(nt,np,t,p)

% gsw_gibbs_ice                     Gibbs energy of ice and its derivatives
% =========================================================================
%
% USAGE:
%  gibbs_ice = gsw_gibbs_ice(nt,np,t,p)
%
% DESCRIPTION:
%  Ice specific Gibbs energy and derivatives up to order 2.
%
% INPUT:
%  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
%  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%   
% OUTPUT:
%  gibbs_ice = Specific Gibbs energy of ice or its derivatives.
%            The Gibbs energy (when nt = np = 0) has units of:     [ J/kg ]
%            The temperature derivatives are output in units of: 
%                                                      [ (J/kg) (K)^(-nt) ]
%            The pressure derivatives are output in units of:
%                                                     [ (J/kg) (Pa)^(-np) ]
%            The mixed derivatives are output in units of:
%                                           [ (J/kg) (K)^(-nt) (Pa)^(-np) ]
%  Note. The derivatives are taken with respect to pressure in Pa, not
%    withstanding that the pressure input into this routine is in dbar.
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IAPWS, 2009: Revised release on the Equation of State 2006 for H2O Ice 
%   Ih. The International Association for the Properties of Water and 
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.  
%
%  Reference page in Help browser
%       <a href="matlab:doc gsw_gibbs_ice">doc gsw_gibbs_ice</a>
%  Note that this reference page includes the code contained in 
%  gsw_gibbs_ice.  We have opted to encode this programme as it is a global
%  standard and such we cannot allow anyone to change it.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

rec_Pt = 1.634903221903779e-3;   % 1./Pt, where Pt = 611.657;  Experimental 
                                 % triple-point pressure in Pa.

Tt = 273.16;  % Triple-point temperature, kelvin (K).  
rec_Tt = 3.660858105139845e-3;   % 1/Tt = 3.660858105139845e-3; 

T = t + 273.15; % The input temperature t is in-situ temperature in
                % units of degrees Celcius.  T is the in-situ Absolute 
                % Temperature of the ice in degrees kelvin (K).  
tau = T.*rec_Tt;

db2Pa = 1e4;
dzi = db2Pa.*p.*rec_Pt;

g00 = -6.32020233335886e5;
g01 =  6.55022213658955e-1;
g02 = -1.89369929326131e-8;
g03 =  3.3974612327105304e-15;
g04 = -5.564648690589909e-22;

s0 = -3.32733756492168e3;

t1 = (3.68017112855051e-2 + 5.10878114959572e-2i);
t2 = (3.37315741065416e-1 + 3.35449415919309e-1i);

r1 = (4.47050716285388e1 + 6.56876847463481e1i);
r20	= (-7.25974574329220e1 - 7.81008427112870e1i);
r21	= (-5.57107698030123e-5 + 4.64578634580806e-5i);
r22	= (	2.34801409215913e-11 - 2.85651142904972e-11i);

if nt == 0 & np == 0
    
    tau_t1 = tau./t1;
    sqtau_t1 = tau_t1.*tau_t1;
    tau_t2 = tau./t2;
    sqtau_t2 = tau_t2.*tau_t2;
    
    g0 = g00 + dzi.*(g01 + dzi.*(g02 + dzi.*(g03 + g04.*dzi)));

    r2 = r20 + dzi.*(r21 + r22.*dzi);
         
    g = r1.*(tau.*log((1 + tau_t1)./(1 - tau_t1)) + t1.*(log(1 - sqtau_t1) - sqtau_t1)) ...
       + r2.*(tau.*log((1 + tau_t2)./(1 - tau_t2)) + t2.*(log(1 - sqtau_t2) - sqtau_t2));
   
    gibbs_ice = g0 - Tt.*(s0.*tau - real(g));
    
elseif nt == 1 & np == 0
    
    tau_t1 = tau./t1;
    tau_t2 = tau./t2;
    
    r2 = r20 + dzi.*(r21 + r22.*dzi);
    
    g = r1.*(log((1 + tau_t1)./(1 - tau_t1)) - 2.*tau_t1) ...
        + r2.*(log((1 + tau_t2)./(1 - tau_t2)) - 2.*tau_t2);
    
    gibbs_ice = -s0 + real(g);
        
elseif nt == 0 & np == 1
    
    tau_t2 = tau./t2;
    sqtau_t2 = tau_t2.*tau_t2;
    
    g0p = rec_Pt.*(g01 + dzi.*(2.*g02 + dzi.*(3.*g03 + 4.*g04.*dzi)));
    
    r2p = rec_Pt.*(r21 + 2.*r22.*dzi);
        
    g = r2p.*(tau.*log((1 + tau_t2)./(1 - tau_t2)) + t2.*(log(1 - sqtau_t2) ...
        - sqtau_t2));
    
    gibbs_ice = g0p + Tt.*real(g);

elseif nt == 1 & np == 1 
    
    tau_t2 = tau./t2;

    r2p = rec_Pt.*(r21 + 2.*r22.*dzi); 
    
    g = r2p.*(log((1 + tau_t2)./(1 - tau_t2)) - 2.*tau_t2);
    
    gibbs_ice = real(g);
    
elseif nt == 2 & np == 0
    
    r2 = r20 + dzi.*(r21 + r22.*dzi);
    
    g = r1.*(1./(t1 - tau) + 1./(t1 + tau) - 2./t1) ...
        + r2.*(1./(t2 - tau) + 1./(t2 + tau) - 2./t2);
    
    gibbs_ice = rec_Tt.*real(g);
    
elseif nt == 0 & np == 2
   
    sqrec_Pt = rec_Pt.*rec_Pt;
    
    tau_t2 = tau./t2;
    sqtau_t2 = tau_t2.*tau_t2;
    
    g0pp = sqrec_Pt.*(2.*g02 + dzi.*(6.*g03 + 12.*g04.*dzi));
    
    r2pp = 2.*r22.*sqrec_Pt;
       
    g = r2pp.*(tau.*log((1 + tau_t2)./(1 - tau_t2)) + t2.*(log(1 - sqtau_t2) ...
       - sqtau_t2));

   gibbs_ice = g0pp + Tt.*real(g);
    
end
    
end