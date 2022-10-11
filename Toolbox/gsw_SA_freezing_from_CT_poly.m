function SA_freezing = gsw_SA_freezing_from_CT_poly(CT,p,saturation_fraction)

% gsw_SA_freezing_from_CT_poly             Absolute Salinity of seawater at
%                                                 the freezing point (poly)
%==========================================================================
%
% USAGE:
%  SA_freezing = gsw_SA_freezing_from_CT_poly(CT,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.  
%  That is, the output is the Absolute Salinity of seawater, with the 
%  fraction saturation_fraction of dissolved air, that is in equilibrium 
%  with ice at Conservative Temperature CT and pressure p.  If the input 
%  values are such that there is no positive value of Absolute Salinity for
%  which seawater is frozen, the output, SA_freezing, is put equal to NaN.
%
% INPUT:
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction  =  the saturation fraction of dissolved air in 
%                          seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where CT is MxN.
%
% OUTPUT:
%  SA_freezing =  Absolute Salinity of seawater when it freezes, for 
%                 given input values of Conservative Temperature
%                 pressure and air saturation fraction.            [ g/kg ]              
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See section 3.33 of this TEOS-10 Manual.  
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_SA_freezing_from_CT_poly:  Requires either two or three inputs')
end 

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_SA_freezing_from_CT_poly: saturation_fraction MUST be between zero and one.')
end

[mt,nt] = size(CT);
[mp,np] = size(p);
[map,nap] = size(saturation_fraction);

if (mp == 1) & (np == 1)                   % p scalar - fill to size of CT
    p = p*ones(size(CT));
elseif (nt == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,mt), :);                          % copy down each column.
elseif (mt == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,nt));                            % copy across each row.
elseif (nt == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,mt), :);                          % copy down each column.
elseif (mt == mp) & (nt == np)
    % ok
else
    error('gsw_SA_freezing_from_CT_poly: Inputs array dimensions arguments do not agree')
end %if

if (map == 1) & (nap == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(CT));         % fill to size of CT
elseif (nt == nap) & (map == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,mt), :);      % copy down each column.
elseif (mt == map) & (nap == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,nt));        % copy across each row.
elseif (nt == map) & (nap == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,mt), :);      % copy down each column.
elseif (mt == map) & (nt == nap)
    % ok
else
    error('gsw_SA_freezing_from_CT_poly: Inputs array dimensions arguments do not agree')
end %if

if mt == 1
    CT = CT.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

aa = 0.014289763856964;           % Note that aa = 0.502500117621/35.16504.
b = 0.057000649899720;

c0  =  0.017947064327968736;
%
c1 =  -6.076099099929818;
c2 =   4.883198653547851;
c3 =  -11.88081601230542;
c4 =   13.34658511480257;
c5 =  -8.722761043208607;
c6 =   2.082038908808201;
%
c7 =  -7.389420998107497;
c8 =  -2.110913185058476;
c9 =   0.2295491578006229; 
%  
c10 = -0.9891538123307282;
c11 = -0.08987150128406496;
c12 =  0.3831132432071728;
c13 =  1.054318231187074;
c14 =  1.065556599652796;
c15 = -0.7997496801694032;
c16 =  0.3850133554097069;
c17 = -2.078616693017569;
c18 =  0.8756340772729538;
c19 = -2.079022768390933;
c20 =  1.596435439942262;
c21 =  0.1338002171109174;
c22 =  1.242891021876471;

p0  =  2.570124672768757e-1;
p1  = -1.917742353032266e+1;
p2  = -1.413382858617969e-2;
p3  = -5.427484830917552e-1;
p4  = -4.126621135193472e-4;
p5  = -4.176407833276121e-7;
p6  =  4.688217641883641e-5;
p7  = -3.039808885885726e-8;
p8  = -4.990118091261456e-11;
p9  = -9.733920711119464e-9;
p10 = -7.723324202726337e-12;
p11 =  7.121854166249257e-16;
p12 =  1.256474634100811e-12;
p13 =  2.105103897918125e-15;
p14 =  8.663811778227171e-19;

p_r = p.*1e-4;

% Form the first estimate of SA_freezing from a polynomial in CT and p_r. 
SA = -(CT + 9*p_r)./0.06;       % A rough estimate to get the saturated CT.

SA(SA < 0) = 0;

% CTsat is the estimated value of CT if the seawater were saturated with
% dissolved air, recognizing that it actually has the air fraction
% saturation_fraction; see McDougall et al. (2014).   
CTsat = CT ...
    - (1-saturation_fraction).*(1e-3).*(2.4-aa.*SA).*(1+b.*(1-SA./35.16504));

SA = p0 + p.*(p2 + p4*CTsat + p.*(p5 + CTsat.*(p7 + p9*CTsat) ...
    + p.*(p8  + CTsat.*(p10 + p12*CTsat) + p.*(p11 + p13*CTsat + p14*p)))) ...
    + CTsat.*(p1 + CTsat.*(p3 + p6*p));

CT_freezing_zero_SA = c0 + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) ...
                         - saturation_fraction.*(2.4e-3).*(1 + b);
                     
% Find CT > CT_freezing_zero_SA.  If this is the case, the input values
% represent seawater that is not frozen (at any positive SA). 
[Itw] = find(CT > CT_freezing_zero_SA);       % Itw stands for "I_too_warm"
 if ~isempty(Itw)
    SA(Itw) = NaN; 
 end

% Find -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA  
% with one based on (CT_freezing_zero_SA - CT).
SA_cut_off = 2.5; % This is the band of SA within +- 2.5 g/kg of SA = 0, 
%                   which we treat differently in calculating the initial
%                   values of both SA and dCT_dSA. 
[Ico] = find(abs(SA) < SA_cut_off);
SA(SA < 0 & SA >= -SA_cut_off) = 0;

% Find SA < -SA_cut_off, set them to NaN.
SA(SA < -SA_cut_off) = NaN;

%--------------------------------------------------------------------------
% Form the first estimate of dCT_dSA, the derivative of CT with respect 
% to SA at fixed p.  
%--------------------------------------------------------------------------
SA_r = 0.01*SA;
x = sqrt(SA_r);
dCT_dSA = (c1 + x.*(1.5*c2  + x.*(2*c3  + x.*(2.5*c4  + x.*(3*c5  + 3.5*c6*x)))) ...
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22*x)))) ... 
       + p_r.*(c12 + x.*(1.5*c14 + x.*(2*c17 + 2.5*c20*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21*x))))).*1e-2...
        - saturation_fraction.*(1e-3).*(-0.018994561378548 - SA.*4.632588654871302e-05);

% Now replace the estimate of SA with the one based on 
% (CT_freezing_zero_SA - CT) when (abs(SA) < SA_cut_off). 
if ~isempty(Ico)
    SA(Ico) = (CT(Ico) - CT_freezing_zero_SA(Ico))./dCT_dSA(Ico);
end

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson method to solve the root of 
% CT_freezing = CT for SA. 
%--------------------------------------------------------------------------
Number_of_Iterations = 2;
for I_iter = 1:Number_of_Iterations

%--------------------------------------------------------------------------
% CT_freezing temperature function evaluation (the forward function 
% evaluation), being the same as gsw_CT_freezing_poly(SA,p,saturation_fraction). 
%--------------------------------------------------------------------------
SA_r = 0.01*SA;
x = sqrt(SA_r);
SA_old = SA;
CT_freezing = c0 + SA_r.*(c1 + x.*(c2 + x.*(c3 + x.*(c4 + x.*(c5 + c6.*x))))) ...
    + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) + SA_r.*p_r.*(c10 + p_r.*(c12 + p_r.*(c15 ...
    + c21.*SA_r)) + SA_r.*(c13 + c17.*p_r + c19.*SA_r) + x.*(c11 + p_r.*(c14 ...
    + c18.*p_r)  + SA_r.*(c16 + c20.*p_r + c22.*SA_r))) ...
    - saturation_fraction.*(1e-3).*(2.4 - aa.*SA).*(1 + b.*(1 - SA./35.16504));

SA = SA_old - (CT_freezing - CT)./dCT_dSA;
                % This is the half-way point of the modified Newton-Raphson 
                % method of McDougall and Wotherspoon (2014). 
SA_r = 0.5*0.01*(SA + SA_old); % This is now the mean value of SA and SA_old. 
x = sqrt(SA_r);
dCT_dSA = (c1 + x.*(1.5*c2  + x.*(2*c3  + x.*(2.5*c4  + x.*(3*c5  + 3.5*c6*x)))) ...
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22*x)))) ... 
       + p_r.*(c12 + x.*(1.5*c14 + x.*(2*c17 + 2.5*c20*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21*x))))).*1e-2...
        - saturation_fraction.*(1e-3).*(-0.018994561378548 - SA.*4.632588654871302e-05);

SA = SA_old - (CT_freezing - CT)./dCT_dSA;
               % This is the end of one full iteration of the modified 
               % Newton-Raphson method of McDougall and Wotherspoon (2013). 
end

%--------------------------------------------------------------------------
% The following lines of code, if implemented, calculates the error of 
% this function in terms of Conservative Temperature, CT_error.  
% With Number_of_Iterations = 1, the maximum error in CT is 2x10^-7 C.
% With Number_of_Iterations = 2, the maximum error in CT is 7x10^-15 C, 
% which is the machine precision of the computer. 
% Number_of_Iterations = 2 is what we recommend. 
%
% SA_r = 0.01*SA;
% x = sqrt(SA_r);
% CT_freezing = c0 ...
%  + SA_r.*(c1 + x.*(c2 + x.*(c3 + x.*(c4 + x.*(c5 + c6.*x))))) ...
%  + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) ...
%  + SA_r.*p_r.*(c10 + p_r.*(c12 + p_r.*(c15 + c21.*SA_r)) + SA_r.*(c13 + c17.*p_r + c19.*SA_r) ...
%  + x.*(c11 + p_r.*(c14 + c18.*p_r)  + SA_r.*(c16 + c20.*p_r + c22.*SA_r))) ...
%  - saturation_fraction.*(1e-3).*(2.4 - aa.*SA).*(1 + b.*(1 - SA./35.16504));
% 
% CT_error = abs(CT_freezing - CT);
% 
% CT_error(p > 10000 | SA > 120 | ...
%     p + SA.*71.428571428571402 > 13571.42857142857) = NaN;
%
%--------------------This is the end of the error calculation--------------

SA_freezing = SA;

%Find SA that are out of range, set them to NaN. 
SA_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if ~isempty(Itw)
    SA_freezing(Itw) = NaN; % If the CT input is too warm, then there is 
end %            no (positive) value of SA that represents frozen seawater. 

if transposed
    SA_freezing = SA_freezing.';
end

end