function enthalpy_SSO_0_CT25 = gsw_enthalpy_SSO_0_CT25(p)
    
% gsw_h_SSO_0_CT25                          enthalpy_CT25 at (SSO,CT = 0,p)
%                                                            (25-term eqn.)
%==========================================================================
%  This function calcualtes enthalpy at the Standard Ocean Salinty, SSO, 
%  and at a Conservative Temperature of zero degrees C, as a function of
%  pressure, p, in dbar, using a streamlined version of the 25-term CT 
%  version of the Gibbs function, that is, a streamlined version of the 
%  code "gsw_enthalpy_CT25(SA,CT,p)".
%==========================================================================

db2Pa = 1e4; 
SSO = 35.16504*ones(size(p));

a0 = 1 + SSO.*(2.0777716085618458e-3 + sqrt(SSO).*3.4688210757917340e-6);
a1 = 6.8314629554123324e-6;

b0 = 9.9984380290708214e2 + SSO.*(2.8925731541277653e0 + ...
    SSO.*1.9457531751183059e-3);

b1 = 0.5*(1.1930681818531748e-2 + SSO.*5.9355685925035653e-6);
b2 = -2.5943389807429039e-8;
b1sq = b1.*b1; 
sqrt_disc = sqrt(b1sq - b0.*b2);
A = b1 - sqrt_disc;
B = b1 + sqrt_disc;

part = (a0.*b2 - a1.*b1)./(b2.*(B - A));

enthalpy_SSO_0_CT25 = db2Pa.*((a1./(2*b2)).*...
    log(1 + p.*(2*b1 + b2.*p)./b0) + ...
    part.*log(1 + (b2.*p.*(B - A))./(A.*(B + b2.*p))));

end
