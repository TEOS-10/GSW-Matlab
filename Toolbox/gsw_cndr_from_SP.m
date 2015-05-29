function R = gsw_cndr_from_SP(SP,t,p)

% gsw_cndr_from_SP               conductivity ratio from Practical Salinity
%==========================================================================
%
% USAGE:  
%  R = gsw_cndr_from_SP(SP,t,p)
%
% DESCRIPTION:
%  Calculates conductivity ratio (R) from (SP,t,p) using PSS-78.  Note that
%  the PSS-78 algorithm for Practical Salinity is only valid in the 
%  range 2 < SP < 42. 
%
% INPUT:
%  SP  =  Practical Salinity  (PSS-78)                         [ unitless ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar) 
%
%  SP & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SP & t are MxN.
%
% OUTPUT:
%  R   =  conductivity ratio                                   [ unitless ]
%
% AUTHOR:  
%  Phil Morgan              [ help_gsw@csiro.au ]
%
% MODIFIED:
%  25th June 1999. Lindsay Pender, Fixed transpose of row vectors.
%  12 December 2003. Lindsay Pender, Converted to ITS-90.
%  16th August 2010. Paul Barker, Regrouped.
%
% VERSION NUMBER: 2.0 (16th August, 2010)
%
% REFERENCES:
%  Fofonoff, P. and R.C .Millard, Jr, 1983: Algorithms for computation of
%   fundamental properties of seawater, 1983. 
%   Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.
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
    error('gsw_cndr_from_SP: Must have 3 input arguments')
end %if

% These few lines ensure that SA is non-negative.
[I_neg_SP] = find(SP < 0);
if ~isempty(I_neg_SP)
    error('gsw_cndr_from_SP: SP must be non-negative!')
end

[ms,ns] = size(SP);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_cndr_from_SP: SP and t must have same dimensions')
end

if mp==1 & np==1                       % p is a scalar.  Fill to size of SP
    p = p(1)*ones(ms,ns);
elseif np==ns & mp==1                   % p is row vector, 
    p = p(ones(1,ms),:);                %   copy down each column.
elseif mp==ms & np==1                   % p is column vector,
    p = p(:,ones(1,ns));                %   copy across each row.
elseif mp==ms & np==ns               
    % ok
else
    error('gsw_cndr_from_SP: p has wrong dimensions')
end %if

if ms == 1
    SP = SP';
    t = t';
    p = p';
    [ms,ns] = size(SP);
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

t68 = t.*1.00024;
d_t68 = t68 - 15;
Rx = nan(size(SP));

for i = 1:ms
    for j = 1:ns
     %---------------------------------------------------------------------
     %DO A NEWTON-RAPHSON ITERATION FOR INVERSE INTERPOLATION OF Rt FROM SP
     %---------------------------------------------------------------------
        SP_ij = SP(i,j);                                   % SP in the loop
        d_t68_ij = d_t68(i,j);                            % t68 in the loop
        Rx_ij = abs(SP_ij/35.0);             % first guess at Rx = sqrt(Rt)
        Rtx = Rx_ij;

        d_S = (d_t68_ij ./ (1+k*d_t68_ij) ) .* ...
            ( b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
        S = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx;
        SP_Inc = S + d_S;
   
        Iloop    = 0;
        end_ij = 0;
        
        while ~end_ij
            Rtx = sqrt(Rx_ij);
            dS =  a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5.*Rtx).*Rtx).*Rtx).*Rtx + ...
                ((d_t68_ij)./(1+k*(d_t68_ij)))* ...
                (b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5.*Rtx).*Rtx).*Rtx).*Rtx);
            Rx_ij = Rx_ij + (SP_ij - SP_Inc)./dS;
                         
            Rtx = Rx_ij;
            d_S = (d_t68_ij ./ (1+k*d_t68_ij) ) .* ...
                (b0 + (b1 + (b2+ (b3 + (b4 + b5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx);
            S = a0 + (a1 + (a2 + (a3 + (a4 + a5.*Rtx).*Rtx).*Rtx).*Rtx).*Rtx;
            SP_Inc = S + d_S;            

            Iloop  = Iloop + 1;
            dels = abs(SP_Inc-SP_ij);            
            if (dels>1.0e-10 & Iloop<100)
                end_ij = 0;
            else
                end_ij = 1;
            end %if
        end %while
        
        Rx(i,j) = Rx_ij;
        
    end %for j
end %for i

%--------------------------------------------------------------------------
% ONCE Rt FOUND, CORRESPONDING TO EACH (SP,t) EVALUATE R. 
% Eqn 4, p.8 (Unesco, 1983)
%--------------------------------------------------------------------------

A  = (d3 + d4.*t68);
B  = 1 + d1.*t68 + d2.*t68.^2;
C  = p.*(e1 + e2.*p + e3.*p.^2);

Rt    = Rx.*Rx;                                % Eqn 6, p.9 (UNESCO, 1983).

rt = c0 + (c1 + (c2 + (c3 + c4.*t68).*t68).*t68).*t68;

D     = B - A.*rt.*Rt;
E     = rt.*Rt.*A.*(B+C);
Ra     = sqrt(abs(D.^2+4*E)) - D;
R     = 0.5*Ra./A;

if transposed
    R = R';
end

end

