function SP = gsw_SP_from_cndr(R,t,p)

% gsw_SP_from_cndr                     Practical Salinity from conductivity
%==========================================================================
%
% USAGE: 
%   SP = gsw_SP_from_cndr(R,t,p)
%
% DESCRIPTION:
%  Calculates Practical Salinity from conductivity ratio (R), using the
%  PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical Salinity
%  is only valid in the range 2 < SP < 42.  The output, SP, of this 
%  function is constrained to be non-negative.    
%  
% INPUT:
%   R     =   Conductivity ratio                               [ unitless ]
%   t     =   in-situ temperature (ITS-90)                        [ deg C ]
%   p     =   sea pressure                                         [ dbar ]
%           ( ie. absolute pressure - 10.1325 dbar )
%
%  t & p may have dimensions 1x1 or Mx1 or 1xN or MxN, where R is MxN.
%
% OUTPUT:
%   SP    =   Practical Salinity on the PSS-78 scale           [ unitless ]
%
% AUTHOR:  
%   17th April 1993.  Phil Morgan          [ help_gsw@csiro.au ]
%
% MODIFIED:
%   12th December 2003. Lindsay Pender, Converted to ITS-90.
%   28th July 2010. by Paul Barker  
%
% VERSION NUMBER: 2.0 (3rd August, 2010)
%
% REFERENCES:
%  Fofonoff, P. and R.C. Millard Jr. 1983: Algorithms for computation of 
%   fundamental properties of seawater. Unesco Tech. Pap. in Mar. Sci., 44,
%   53 pp.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix E of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_SP_from_cndr.m:  Requires three input arguments')
end %if

[mc,nc] = size(R);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt == 1) & (nt == 1)              % t scalar - fill to size of R
    t = t*ones(size(R));
elseif (nc == nt) & (mt == 1)         % t is row vector,
    t = t(ones(1,mc), :);              % copy down each column.
elseif (mc == mt) & (nt == 1)         % t is column vector,
    t = t(:,ones(1,nc));               % copy across each row.
elseif (mc == mt) & (nc == nt)
    % ok
else
    error('gsw_SP_from_cndr.m: Inputs array dimensions arguments do not agree')
end %if

if (mp == 1) & (np == 1)              % p scalar - fill to size of R
    p = p*ones(size(R));
elseif (nc == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,mc), :);              % copy down each column.
elseif (mc == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nc));               % copy across each row.
elseif (mc == mp) & (nc == np)
    % ok
else
    error('gsw_SP_from_cndr.m: Inputs array dimensions arguments do not agree')
end %if

if mc == 1
    R = R';
    t = t';
    p = p';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

t68 = t * 1.00024;
d_t68 = t68 - 15;

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

k  =  0.0162;

rt = c0 + (c1 + (c2 + (c3 + c4.*t68).*t68).*t68).*t68;
Rp = 1 + ( p.*(e1 + e2.*p + e3.*p.*p) ) ./ ...
      (1 + d1.*t68 + d2.*t68.*t68 + (d3 + d4.*t68).*R);
Rt = R./(Rp.*rt);
Rtx   = sqrt(Rt);

d_S = (d_t68 ./ (1 + k*d_t68)) .* ...
        (b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
S = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx;
SP = S + d_S;
% 
% These few lines ensure that SP is non-negative. 
[I_neg_SP] = find(SP < 0);
if ~isempty(I_neg_SP)
    SP(I_neg_SP) = 0;
end

if transposed
    SP = SP';
end

end
