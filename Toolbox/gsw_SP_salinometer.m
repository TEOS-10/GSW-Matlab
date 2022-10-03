function SP = gsw_SP_salinometer(Rt,t)

% gsw_SP_salinometer       Practical Salinity from a laboratory salinometer
%==========================================================================
%
% USAGE: 
%  SP = gsw_SP_salinometer(Rt,t)
%
% DESCRIPTION:
%  Calculates Practical Salinity SP from a salinometer, primarily using the 
%  PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity
%  is only valid in the range 2 < SP < 42.  If the PSS-78 algorithm 
%  produces a Practical Salinity that is less than 2 then the Practical 
%  Salinity is recalculated with a modified form of the Hill et al. (1986) 
%  formula.  The modification of the Hill et al. (1986) expression is to 
%  ensure that it is exactly consistent with PSS-78 at SP = 2. 
%
%  A laboratory salinometer has the ratio of conductivities, Rt, as an 
%  output, and the present function uses this conductivity ratio and the
%  temperature t of the salinometer bath as the two input variables. 
%  
% INPUT:
%  Rt  =  C(SP,t_68,0)/C(SP=35,t_68,0)                         [ unitless ]
%  t   =  temperature of the bath of the salinometer,
%         measured on the ITS-90 scale (ITS-90)                   [ deg C ]
%
% OUTPUT:
%  SP  =  Practical Salinity on the PSS-78 scale               [ unitless ]
%
%  t may have dimensions 1x1 or Mx1 or 1xN or MxN, where Rt is MxN.
%
% AUTHOR:  
%  Paul Barker, Trevor McDougall and Rich Pawlowicz    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation of 
%   fundamental properties of seawater. Unesco Tech. Pap. in Mar. Sci., 44,
%   53 pp.
%
%  Hill, K.D., T.M. Dauphinee & D.J. Woods, 1986: The extension of the 
%   Practical Salinity Scale 1978 to low salinities. IEEE J. Oceanic Eng.,
%   11, 109 - 112.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See appendix E of this TEOS-10 Manual, and in particular, 
%     Eqns. (E.2.1) and (E.2.6). 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_SP_salinometer:  Requires two input arguments')
end %if

[mc,nc] = size(Rt);
[mt,nt] = size(t);

if (mt == 1) & (nt == 1)               % t scalar - fill to size of Rt
    t = t*ones(size(Rt));
elseif (nc == nt) & (mt == 1)          % t is row vector,
    t = t(ones(1,mc), :);               % copy down each column.
elseif (mc == mt) & (nt == 1)          % t is column vector,
    t = t(:,ones(1,nc));                % copy across each row.
elseif (nc == mt) & (np == 1)          % t is a transposed row vector,
    t = t.';                                         % transposed then
    t = t(ones(1,mc), :);                    % copy down each column.
elseif (mc == mt) & (nc == nt)
    % ok
else
    error('gsw_SP_salinometer: Inputs array dimensions arguments do not agree')
end %if

if mc == 1
    Rt = Rt.';
    t = t.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
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

k  =  0.0162;

t68 = t.*1.00024;
ft68 = (t68 - 15)./(1 + k*(t68 - 15));
 
Rt(Rt < 0) = NaN;

Rtx = sqrt(Rt);

SP = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx + ...
    ft68.*(b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);

% The following section of the code is designed for SP < 2 based on the
% Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
% exactly equal to the PSS-78 algorithm at SP = 2.

if any(SP < 2)
    [I2] = find(SP < 2);
    Hill_ratio = gsw_Hill_ratio_at_SP2(t(I2));  
    x = 400*Rt(I2);
    sqrty = 10*Rtx(I2);   
    part1 = 1 + x.*(1.5 + x) ;
    part2 = 1 + sqrty.*(1 + sqrty.*(1 + sqrty));
    SP_Hill_raw = SP(I2) - a0./part1 - b0.*ft68(I2)./part2;
    SP(I2) = Hill_ratio.*SP_Hill_raw;
end

% This line ensures that SP is non-negative.
SP(SP < 0) = 0;

if transposed
    SP = SP.';
end

end
