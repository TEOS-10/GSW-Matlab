function SA_freezing = gsw_SA_freezing_from_t_poly(t,p,saturation_fraction)

% gsw_SA_freezing_from_t_poly          Absolute Salinity of seawater at the
%                                                     freezing point (poly)
%==========================================================================
%
% USAGE:
%  SA_freezing = gsw_SA_freezing_from_t_poly(t,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.  
%  That is, the output is the Absolute Salinity of seawater, with the 
%  fraction saturation_fraction of dissolved air, that is in equilibrium 
%  with ice at in-situ temperature t and pressure p.  If the input values 
%  are such that there is no positive value of Absolute Salinity for which 
%  seawater is frozen, the output, SA_freezing, is put equal to NaN.  
%
% INPUT:
%  t  =  in-situ Temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where t is MxN.
%
% OUTPUT:
%  SA_freezing  =  Absolute Salinity of seawater when it freezes, for 
%                  given input values of in situ temperature, pressure and 
%                  air saturation fraction.                        [ g/kg ]              
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
   error('gsw_SA_freezing_from_t_poly:  Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_SA_freezing_from_t_poly: saturation_fraction MUST be between zero and one.')
end
   
[mt,nt] = size(t);
[mp,np] = size(p);
[map,nap] = size(saturation_fraction);

if (mp == 1) & (np == 1)                   % p scalar - fill to size of t
    p = p*ones(size(t));
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
    error('gsw_SA_freezing_from_t_poly: Inputs array dimensions arguments do not agree')
end %if

if (map == 1) & (nap == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(t));         % fill to size of t
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
    error('gsw_SA_freezing_from_t_poly: Inputs array dimensions arguments do not agree')
end %if

if mt == 1
    t = t.';
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

p0  =  2.570124672768757e-1;
p1  = -1.917742353032266e1;
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

% Coefficients from gsw_t_freezing_poly
%c0 = 0.002519;
c1 = -5.946302841607319;
c2 =  4.136051661346983;
c3 = -1.115150523403847e1;
c4 =  1.476878746184548e1;
c5 = -1.088873263630961e1;
c6 =  2.961018839640730;
% c7 = -7.433320943962606;
% c8 = -1.561578562479883;
% c9 =  4.073774363480365e-2;
c10 =  1.158414435887717e-2;
c11 = -4.122639292422863e-1;
c12 = -1.123186915628260e-1;
c13 =  5.715012685553502e-1;
c14 =  2.021682115652684e-1;
c15 =  4.140574258089767e-2;
c16 = -6.034228641903586e-1;
c17 = -1.205825928146808e-2;
c18 = -2.812172968619369e-1;
c19 =  1.877244474023750e-2;
c20 = -1.204395563789007e-1;
c21 =  2.349147739749606e-1;
c22 =  2.748444541144219e-3;

p_r = p.*1e-4;

%--------------------------------------------------------------------------
% Form the first estimate of SA_freezing, called SA here, from a polynomial 
% in CT and p_r. 
%--------------------------------------------------------------------------
SA = -(t + 9*p_r)./0.06; % A rough estimate to get the saturated CT.

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

CT = gsw_CT_from_t(SA,t,p);

% CTsat is the estimated value of CT if the seawater were saturated with
% dissolved air, recognizing that it actually has the air fraction
% saturation_fraction; see McDougall et al. (2014).  
CTsat = CT ...
    - (1-saturation_fraction).*(1e-3).*(2.4-aa.*SA).*(1+b.*(1-SA./35.16504));

SA = p0 + p.*(p2 + p4*CTsat + p.*(p5 + CTsat.*(p7 + p9*CTsat) ...
    + p.*(p8  + CTsat.*(p10 + p12*CTsat) + p.*(p11 + p13*CTsat + p14*p)))) ...
    + CTsat.*(p1 + CTsat.*(p3 + p6*p));

t_freezing_zero_SA = gsw_t_freezing_poly(zeros(size(t)),p,saturation_fraction);

% Find t > t_freezing_zero_SA.  If this is the case, the input values
% represent seawater that is not frozen (at any positive SA). 
[Itw] = find(t > t_freezing_zero_SA);         % Itw stands for "I_too_warm"
 if ~isempty(Itw)
    SA(Itw) = NaN; 
 end
 
% Find -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA  
% with one based on (t_freezing_zero_SA - t).
SA_cut_off = 2.5;   % This is the band of SA within +- 2.5 g/kg of SA = 0, 
%                     which we treat differently in calculating the initial
%                     values of both SA and dCT_dSA. 
[Ico] = find(abs(SA) < SA_cut_off);

SA(SA < 0 & SA >= -SA_cut_off) = 0;

% Find SA < -SA_cut_off, set them to NaN.
SA(SA < -SA_cut_off) = NaN;

%---------------------------------------------------------------------------
% Form the first estimate of dt_dSA, the derivative of t with respect 
% to SA at fixed p, using the coefficients, c0 ... c22 from gsw_t_freezing_poly. 
%--------------------------------------------------------------------------
SA_r = 0.01*SA;
x = sqrt(SA_r);
dt_dSA = (c1 + x.*(1.5*c2 + x.*(2*c3 + x.*(2.5*c4 + x.*(3*c5 + 3.5*c6.*x))))  ... 
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22.*x)))) ...
       + p_r.*(c12*+ x.*(1.5*c14 + x.*(2*c17 + 2.5*c20.*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21.*x))))).*1e-2 ...
       + 1.421866717626370e-05.*saturation_fraction;

% Now replace the estimate of SA with the one based on 
% (t_freezing_zero_SA - t) when (abs(SA) < SA_cut_off). 
if ~isempty(Ico)
    SA(Ico) = (t(Ico) - t_freezing_zero_SA(Ico))./dt_dSA(Ico);
end

%---------------------------------------------------------------------------
% Begin the modified Newton-Raphson method to find the root of 
% t_freezing = t for SA. 
%---------------------------------------------------------------------------
Number_of_Iterations = 5;
for I_iter = 1:Number_of_Iterations
    
    SA_old = SA;
    
    t_freezing = gsw_t_freezing_poly(SA_old,p,saturation_fraction);
    
    SA = SA_old - (t_freezing - t)./dt_dSA;
                % This is the half-way point of the modified Newton-Raphson 
                % method of McDougall and Wotherspoon (2014). 
    SA_r = 0.5*0.01*(SA + SA_old); % This is now the mean value of SA and SA_old.
    x = sqrt(SA_r);
    dt_dSA = (c1 + x.*(1.5*c2 + x.*(2*c3 + x.*(2.5*c4 + x.*(3*c5 + 3.5*c6.*x))))  ... 
       + p_r.*(c10 + x.*(1.5*c11 + x.*(2*c13 + x.*(2.5*c16 + x.*(3*c19 + 3.5*c22.*x)))) ...
       + p_r.*(c12*+ x.*(1.5*c14 + x.*(2*c17 + 2.5*c20.*x)) ...
       + p_r.*(c15 + x.*(1.5*c18 + 2*c21.*x))))).*1e-2 ...
       + 1.421866717626370e-05.*saturation_fraction;
    SA = SA_old - (t_freezing - t)./dt_dSA;
               % This is the end of one full iteration of the modified 
               % Newton-Raphson method of McDougall and Wotherspoon (2014). 
end

%--------------------------------------------------------------------------
% The following lines of code, if implemented, calculate the error of 
% this function in terms of in-situ temperature, t.  
% With Number_of_Iterations = 4, the maximum error in t is 3x10^-13 C.
% With Number_of_Iterations = 5, the maximum error in t is 2x10^-14 C, 
% which is the machine precision of the computer. 
% Number_of_Iterations = 5 is what we recommend. 
%
% SA(SA < 0) = 0;  
% 
% t_freezing = gsw_t_freezing_poly(SA,p,saturation_fraction);
%  
% t_error = abs(t_freezing - t);
%  
% t_error(p > 10000 | SA > 120 | ...
%      p + SA.*71.428571428571402 > 13571.42857142857) = NaN;
% 
%-----------------This is the end of the error calculation-----------------

SA_freezing = SA;

%Find SA that are out of range, set them to NaN. 
SA_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if ~isempty(Itw)
    SA_freezing(Itw) = Nan;       % If the t input is too warm, then there is 
%                no (positive) value of SA that represents frozen seawater. 
end

if transposed
    SA_freezing = SA_freezing.';
end

end