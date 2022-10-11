function R = gsw_R_from_SP(SP,t,p)

% gsw_R_from_SP                                  conductivity ratio from SP
%==========================================================================
%
% USAGE:
%  R = gsw_R_from_SP(SP,t,p)
%
% DESCRIPTION:
%  Calculates conductivity ratio from (SP,t,p) using PSS-78 in the range 
%  2 < SP < 42.  If the input Practical Salinity is less than 2 then a 
%  modified form of the Hill et al. (1986) fomula is used for Practical 
%  Salinity.  The modification of the Hill et al. (1986) expression is to
%  ensure that it is exactly consistent with PSS-78 at SP = 2.
%
%  The conductivity ratio returned by this function is consistent with the
%  input value of Practical Salinity, SP, to 2x10^-14 psu over the full 
%  range of input parameters (from pure fresh water up to SP = 42 psu).  
%  This error of 2x10^-14 psu is machine precision at typical seawater 
%  salinities.  This accuracy is achieved by having four different 
%  polynomials for the starting value of Rtx (the square root of Rt) in 
%  four different ranges of SP, and by using one and a half iterations of 
%  a computationally efficient modified Newton-Raphson technique 
%  (McDougall and Wotherspoon, 2013) to find the root of the equation.  
%
%  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
%  Salinity in terms of the conductivity ratio, R, without actually
%  specifying the value of C(35,15,0) (which we currently take to be
%  42.9140 mS cm^-1 (Culkin and Smith, 1980)).
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SP & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SP & t are MxN.
%
% OUTPUT:
%  R  =  conductivity ratio                                    [ unitless ]
%
% AUTHOR:
%  Trevor McDougall, Paul Barker and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Culkin and Smith, 1980:  Determination of the Concentration of Potassium  
%   Chloride Solution Having the Same Electrical Conductivity, at 15C and 
%   Infinite Frequency, as Standard Seawater of Salinity 35.0000 
%   (Chlorinity 19.37394), IEEE J. Oceanic Eng, 5, 22-23.
%
%  Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the 
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   11, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix E of this TEOS-10 Manual. 
%
%  McDougall T.J., and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  Unesco, 1983: Algorithms for computation of fundamental properties of 
%   seawater. Unesco Technical Papers in Marine Science, 44, 53 pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_R_from_SP: Must have 3 input arguments')
end %if

% These few lines ensure that SP is non-negative.
if any(SP < 0)
    error('gsw_R_from_SP: SP must be non-negative!')
end

[ms,ns] = size(SP);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_R_from_SP: SP and t must have same dimensions')
end

if mp==1 & np==1                       % p is a scalar.  Fill to size of SP
    p = p(1)*ones(ms,ns);
elseif np==ns & mp==1                       % p is row vector,
    p = p(ones(1,ms),:);                    %   copy down each column.
elseif mp==ms & np==1                       % p is column vector,
    p = p(:,ones(1,ns));                    %   copy across each row.
elseif (np == ms) & (np == 1)               % p is a transposed row vector,
    p = p.';                                 % transposed then
    p = p(ones(1,ms), :);                   % copy down each column.
elseif mp==ms & np==ns
    % ok
else
    error('gsw_R_from_SP: p has wrong dimensions')
end %if

if ms == 1
    SP = SP.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Setting up the constants
%--------------------------------------------------------------------------

a0 =  0.0080;
a1 = -0.1692;
a2 = 25.3851;
a3 = 14.0941;
a4 = -7.0261;
a5 =  2.7081;

b0 =  0.0005;
b1 = -0.0056;
b2 = -0.0066;
b3 = -0.0375;
b4 =  0.0636;
b5 = -0.0144;

c0 =  0.6766097;
c1 =  2.00564e-2;
c2 =  1.104259e-4;
c3 = -6.9698e-7;
c4 =  1.0031e-9;

d1 =  3.426e-2;
d2 =  4.464e-4;
d3 =  4.215e-1;
d4 = -3.107e-3;

e1 =  2.070e-5;
e2 = -6.370e-10;
e3 =  3.989e-15;

p0 =   4.577801212923119e-3;
p1 =   1.924049429136640e-1;
p2 =   2.183871685127932e-5;
p3 =  -7.292156330457999e-3;
p4 =   1.568129536470258e-4;
p5 =  -1.478995271680869e-6;
p6 =   9.086442524716395e-4;
p7 =  -1.949560839540487e-5;
p8 =  -3.223058111118377e-6;
p9 =   1.175871639741131e-7;
p10 = -7.522895856600089e-5;
p11 = -2.254458513439107e-6;
p12 =  6.179992190192848e-7;
p13 =  1.005054226996868e-8;
p14 = -1.923745566122602e-9;
p15 =  2.259550611212616e-6;
p16 =  1.631749165091437e-7;
p17 = -5.931857989915256e-9;
p18 = -4.693392029005252e-9;
p19 =  2.571854839274148e-10;
p20 =  4.198786822861038e-12;
    
q0 =   5.540896868127855e-5;
q1 =   2.015419291097848e-1;
q2 =  -1.445310045430192e-5;
q3 =  -1.567047628411722e-2;
q4 =   2.464756294660119e-4;
q5 =  -2.575458304732166e-7;
q6 =   5.071449842454419e-3;
q7 =  -9.081985795339206e-5;
q8 =  -3.635420818812898e-6;
q9 =   2.249490528450555e-8;
q10 = -1.143810377431888e-3;
q11 =  2.066112484281530e-5;
q12 =  7.482907137737503e-7;
q13 =  4.019321577844724e-8;
q14 = -5.755568141370501e-10;
q15 =  1.120748754429459e-4;
q16 = -2.420274029674485e-6;
q17 = -4.774829347564670e-8;
q18 = -4.279037686797859e-9;
q19 = -2.045829202713288e-10;
q20 =  5.025109163112005e-12;

r0 =   3.432285006604888e-3;
r1 =   1.672940491817403e-1;
r2 =   2.640304401023995e-5;
r3 =   1.082267090441036e-1;
r4 =  -6.296778883666940e-5;
r5 =  -4.542775152303671e-7;
r6 =  -1.859711038699727e-1;
r7 =   7.659006320303959e-4;
r8 =  -4.794661268817618e-7;
r9 =   8.093368602891911e-9;
r10 =  1.001140606840692e-1;
r11 = -1.038712945546608e-3;
r12 = -6.227915160991074e-6;
r13 =  2.798564479737090e-8;
r14 = -1.343623657549961e-10;
r15 =  1.024345179842964e-2;
r16 =  4.981135430579384e-4;
r17 =  4.466087528793912e-6;
r18 =  1.960872795577774e-8;
r19 = -2.723159418888634e-10;
r20 =  1.122200786423241e-12;

u0 =    5.180529787390576e-3;
u1 =    1.052097167201052e-3;
u2 =    3.666193708310848e-5;
u3 =    7.112223828976632;
u4 =   -3.631366777096209e-4;
u5 =   -7.336295318742821e-7;
u6 =   -1.576886793288888e+2;
u7 =   -1.840239113483083e-3;
u8 =    8.624279120240952e-6;
u9 =    1.233529799729501e-8;
u10 =   1.826482800939545e+3;
u11 =   1.633903983457674e-1;
u12 =  -9.201096427222349e-5;
u13 =  -9.187900959754842e-8;
u14 =  -1.442010369809705e-10;
u15 =  -8.542357182595853e+3;
u16 =  -1.408635241899082;
u17 =   1.660164829963661e-4;
u18 =   6.797409608973845e-7;
u19 =   3.345074990451475e-10;
u20 =   8.285687652694768e-13;

k  =  0.0162;

t68 = t.*1.00024;
ft68 = (t68 - 15)./(1 + k.*(t68 - 15));

x = sqrt(SP);
Rtx = nan(size(SP));

%--------------------------------------------------------------------------
% Finding the starting value of Rtx, the square root of Rt, using four 
% different polynomials of SP and t68.  
%--------------------------------------------------------------------------

if any(SP >= 9)
    [I] = find(SP >= 9);
    Rtx(I) =  p0 + x(I).*(p1 + p4*t68(I) + x(I).*(p3 + p7*t68(I) + x(I).*(p6 ...
        + p11*t68(I) + x(I).*(p10 + p16*t68(I)+ x(I).*p15))))...
        + t68(I).*(p2+ t68(I).*(p5 + x(I).*x(I).*(p12 + x(I).*p17) + p8*x(I) ...
        + t68(I).*(p9 + x(I).*(p13 + x(I).*p18)+ t68(I).*(p14 + p19*x(I) + p20*t68(I)))));
end

if any(SP >= 0.25 & SP < 9)
    [I] = find(SP >= 0.25 & SP < 9);
    Rtx(I) =  q0 + x(I).*(q1 + q4*t68(I) + x(I).*(q3 + q7*t68(I) + x(I).*(q6 ...
        + q11*t68(I) + x(I).*(q10 + q16*t68(I)+ x(I).*q15))))...
        + t68(I).*(q2+ t68(I).*(q5 + x(I).*x(I).*(q12 + x(I).*q17) + q8*x(I) ...
        + t68(I).*(q9 + x(I).*(q13 + x(I).*q18)+ t68(I).*(q14 + q19*x(I) + q20*t68(I)))));
end

if any(SP >= 0.003 & SP < 0.25)
    [I] = find(SP >= 0.003 & SP < 0.25);
    Rtx(I) =  r0 + x(I).*(r1 + r4*t68(I) + x(I).*(r3 + r7*t68(I) + x(I).*(r6 ...
        + r11*t68(I) + x(I).*(r10 + r16*t68(I)+ x(I).*r15))))...
        + t68(I).*(r2+ t68(I).*(r5 + x(I).*x(I).*(r12 + x(I).*r17) + r8*x(I) ...
        + t68(I).*(r9 + x(I).*(r13 + x(I).*r18)+ t68(I).*(r14 + r19*x(I) + r20*t68(I)))));
end

if any(SP < 0.003)
    [I] = find(SP < 0.003);
    Rtx(I) =  u0 + x(I).*(u1 + u4*t68(I) + x(I).*(u3 + u7*t68(I) + x(I).*(u6 ...
        + u11*t68(I) + x(I).*(u10 + u16*t68(I)+ x(I).*u15))))...
        + t68(I).*(u2+ t68(I).*(u5 + x(I).*x(I).*(u12 + x(I).*u17) + u8*x(I) ...
        + t68(I).*(u9 + x(I).*(u13 + x(I).*u18)+ t68(I).*(u14 + u19*x(I) + u20*t68(I)))));
end

%--------------------------------------------------------------------------
% Finding the starting value of dSP_dRtx, the derivative of SP with respect
% to Rtx.  
%--------------------------------------------------------------------------
dSP_dRtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5.*Rtx).*Rtx).*Rtx).*Rtx ...
    + ft68.*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5.*Rtx).*Rtx).*Rtx).*Rtx);
if any(SP < 2)
    [I2] = find(SP < 2);
    x = 400.*(Rtx(I2).*Rtx(I2));
    sqrty = 10.*Rtx(I2);
    part1 = 1 + x.*(1.5 + x) ;
    part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
    Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));
    dSP_dRtx(I2) = dSP_dRtx(I2)...
        + a0.*800.*Rtx(I2).*(1.5 + 2*x)./(part1.*part1)...
        + b0.*ft68(I2).*(10 + sqrty.*(20 + 30.*sqrty))./(part2.*part2);
    dSP_dRtx(I2) = Hill_ratio.*dSP_dRtx(I2);
end

%--------------------------------------------------------------------------
% One iteration through the modified Newton-Raphson method achieves an
% error in Practical Salinity of about 10^-12 for all combinations of the
% inputs.  One and a half iterations of the modified Newton-Raphson method  
% achevies a maximum error in terms of Practical Salinity of better than 
% 2x10^-14 everywhere. 
%
% We recommend one and a half iterations of the modified Newton-Raphson
% method. 
%
% Begin the modified Newton-Raphson method.  
%--------------------------------------------------------------------------

    SP_est = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx ...
        + ft68.*(b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
    if any(SP_est < 2)
        [I2] = find(SP_est < 2);
        x = 400.*(Rtx(I2).*Rtx(I2));
        sqrty = 10.*Rtx(I2);
        part1 = 1 + x.*(1.5 + x) ;
        part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
        SP_Hill_raw = SP_est(I2) - a0./part1 - b0.*ft68(I2)./part2;
        Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));
        SP_est(I2) = Hill_ratio.*SP_Hill_raw;
    end
 
    Rtx_old = Rtx;
    Rtx = Rtx_old - (SP_est - SP)./dSP_dRtx;
    Rtxm = 0.5*(Rtx + Rtx_old);      % This mean value of Rtx, Rtxm, is the  
%               value of Rtx at which the derivative dSP_dRtx is evaluated.
    
    dSP_dRtx =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5.*Rtxm).*Rtxm).*Rtxm).*Rtxm ...
        + ft68.*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5.*Rtxm).*Rtxm).*Rtxm).*Rtxm);
    
    if any(SP_est < 2)
        [I2] = find(SP_est < 2);
        x = 400.*(Rtxm(I2).*Rtxm(I2));
        sqrty = 10.*Rtxm(I2);
        part1 = 1 + x.*(1.5 + x) ;
        part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
        dSP_dRtx(I2) = dSP_dRtx(I2)...
            + a0.*800.*Rtxm(I2).*(1.5 + 2*x)./(part1.*part1)...
            + b0.*ft68(I2).*(10 + sqrty.*(20 + 30.*sqrty))./(part2.*part2);
        Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));
        dSP_dRtx(I2) = Hill_ratio.*dSP_dRtx(I2);
    end

%--------------------------------------------------------------------------
% The line below is where Rtx is updated at the end of the one full 
% iteration of the modified Newton-Raphson technique.
%--------------------------------------------------------------------------
    Rtx = Rtx_old - (SP_est - SP)./dSP_dRtx;
%--------------------------------------------------------------------------
% Now we do another half iteration of the modified Newton-Raphson  
% technique, making a total of one and a half modified N-R iterations.
%-------------------------------------------------------------------------- 
SP_est = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx ...
        + ft68.*(b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
    if any(SP_est < 2)
        [I2] = find(SP_est < 2);
        x = 400.*(Rtx(I2).*Rtx(I2));
        sqrty = 10.*Rtx(I2);
        part1 = 1 + x.*(1.5 + x) ;
        part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
        SP_Hill_raw = SP_est(I2) - a0./part1 - b0.*ft68(I2)./part2;
        Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));
        SP_est(I2) = Hill_ratio.*SP_Hill_raw;
    end
    Rtx = Rtx - (SP_est - SP)./dSP_dRtx;
    
%--------------------------------------------------------------------------
% The following lines of code are commented out, but when activated, return
% the error, SP_error, in Rtx (in terms of psu). 
% SP_est = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx ...
%     + ft68.*(b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
% if any(SP_est < 2)
%     [I2] = find(SP_est < 2);
%     x = 400.*(Rtx(I2).*Rtx(I2));
%     sqrty = 10.*Rtx(I2);
%     part1 = 1 + x.*(1.5 + x) ;
%     part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
%     SP_Hill_raw = SP_est(I2) - a0./part1 - b0.*ft68(I2)./part2;
%     Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));
%     SP_est(I2) = Hill_ratio.*SP_Hill_raw;
% end
% 
% SP_error = abs(SP - SP_est);
%
%--------------This is the end of the error testing------------------------


%--------------------------------------------------------------------------
% Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
%--------------------------------------------------------------------------
Rt = Rtx.*Rtx;
A  = d3 + d4.*t68;
B  = 1 + d1.*t68 + d2.*t68.^2;
C  = p.*(e1 + e2.*p + e3.*p.^2); 
% rt_lc (i.e. rt_lower_case) corresponds to rt as defined in 
% the UNESCO 44 (1983) routines.
rt_lc = c0 + (c1 + (c2 + (c3 + c4.*t68).*t68).*t68).*t68;

D  = B - A.*rt_lc.*Rt;
E  = rt_lc.*Rt.*A.*(B + C);
Ra = sqrt(D.^2 + 4*E) - D;
R  = 0.5*Ra./A;

if transposed
    R = R';
end

end
