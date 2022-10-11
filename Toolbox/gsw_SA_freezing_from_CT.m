function SA_freezing = gsw_SA_freezing_from_CT(CT,p,saturation_fraction)

% gsw_SA_freezing_from_CT              Absolute Salinity of seawater at the
%                                                            freezing point
%==========================================================================
%
% USAGE:
%  SA_freezing = gsw_SA_freezing_from_CT(CT,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of seawater at the freezing temperature.  
%  That is, the output is the Absolute Salinity of seawater, with 
%  Conservative Temperature CT, pressure p and the fraction 
%  saturation_fraction of dissolved air, that is in equilibrium 
%  with ice at the same in situ temperature and pressure.  If the input 
%  values are such that there is no positive value of Absolute Salinity for
%  which seawater is frozen, the output, SA_freezing, is made a NaN.
%
% INPUT:
%  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
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
%                 given input values of its Conservative Temperature,
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
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of 
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
   error('gsw_SA_freezing_from_CT:  Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_SA_freezing_from_CT: saturation_fraction MUST be between zero and one.')
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
    error('gsw_SA_freezing_from_CT: Inputs array dimensions arguments do not agree')
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
    error('gsw_SA_freezing_from_CT: Inputs array dimensions arguments do not agree')
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

SA = -(CT + 9e-4*p)./0.06;  % Firstly, this is a very rough estimate of SA 
                            % simply to get the saturated CT, CTsat.

SA(SA < 0) = 0; % This line ensures that SA is not negative

% CTsat is the estimated value of CT if the seawater were saturated with
% dissolved air, recognizing that it actually has the air fraction
% saturation_fraction; see McDougall et al, 2014).  
CTsat = CT ...
    - (1-saturation_fraction).*(1e-3).*(2.4-aa.*SA).*(1+b.*(1-SA./35.16504));

% This is the inital guess of SA using a purpose-built polynomial in CTsat and p.   
SA = p0 + p.*(p2 + p4*CTsat + p.*(p5 + CTsat.*(p7 + p9*CTsat) ...
    + p.*(p8  + CTsat.*(p10 + p12*CTsat) + p.*(p11 + p13*CTsat + p14*p)))) ...
    + CTsat.*(p1 + CTsat.*(p3 + p6*p));

CT_freezing_zero_SA =  gsw_CT_freezing(zeros(size(SA)),p,saturation_fraction);                     
% Find CT > CT_freezing_zero_SA.  If this is the case, the input values
% represent seawater that is not frozen (for any positive SA). 
[Itw] = find(CT > CT_freezing_zero_SA);       % Itw stands for "I_too_warm"
 if ~isempty(Itw)
    SA(Itw) = NaN; 
    p(Itw) = NaN;
 end

% Find -SA_cut_off < SA < SA_cut_off, and replace the above estimate of SA  
% with one based on (CT_freezing_zero_SA - CT).
SA_cut_off = 2.5; % This is the band of SA within +- 2.5 g/kg of SA = 0, 
%                   which we treat differently in calculating the initial
%                   values of SA. 
[Ico] = find(abs(SA) < SA_cut_off);
SA(SA < 0 & SA >= -SA_cut_off) = 0;

% Find SA < -SA_cut_off, set them to NaN.
SA(SA < -SA_cut_off) = NaN;

%--------------------------------------------------------------------------
% Form the first estimate of CTfreezing_SA, the derivative of CT_freezing 
% with respect to SA at fixed p.  
%--------------------------------------------------------------------------
[CTfreezing_SA, dummy] = gsw_CT_freezing_first_derivatives(SA,p,saturation_fraction);
% Now replace the estimate of SA with the one based on 
% (CT_freezing_zero_SA - CT) when (abs(SA) < SA_cut_off). 
if ~isempty(Ico)
    SA(Ico) = (CT(Ico) - CT_freezing_zero_SA(Ico))./CTfreezing_SA(Ico);
end

%--------------------------------------------------------------------------
% Begin the modified Newton-Raphson method to solve  
% f = (CT_freezing - CT) = 0 for SA. 
%--------------------------------------------------------------------------
Number_of_Iterations = 3;
for I_iter = 1:Number_of_Iterations

  SA_old = SA;
  f = gsw_CT_freezing(SA,p,saturation_fraction) - CT;
  SA = SA_old - f./CTfreezing_SA;
               % This is the half-way point of the modified Newton-Raphson   
               % method of McDougall and Wotherspoon (2014). 
  SA_mean = 0.5*(SA + SA_old); 
  [CTfreezing_SA, dummy] = gsw_CT_freezing_first_derivatives(SA_mean,p,saturation_fraction);
  SA = SA_old - f./CTfreezing_SA;
               % This is the end of one full iteration of the modified    
               % Newton-Raphson method of McDougall and Wotherspoon (2014). 
end

%--------------------------------------------------------------------------
% The following lines of code, if implemented, calculates the error of 
% this function in terms of Conservative Temperature, CT_error.  
% With Number_of_Iterations = 3, the maximum error in CT is 3.5x10^-13 C
% and in SA it is 6x10^-12 g/kg, which is the machine precision of the 
% computer. 
%
% CT_freezing = gsw_CT_freezing(SA,p,saturation_fraction);
% 
% CT_error = abs(CT_freezing - CT);
% 
% CT_error(p > 10000 | SA > 120 | ...
%     p + SA.*71.428571428571402 > 13571.42857142857) = NaN;
% if ~isempty(Itw)
%     CT_error(Itw) = NaN;    % If the CT input is too warm, then there is
% end %            no (positive) value of SA that represents frozen seawater.
%
%--------------------This is the end of the error calculation--------------

SA_freezing = SA;

%Find SA that are out of range, set them to NaN. 
SA_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if ~isempty(Itw)
    SA_freezing(Itw) = NaN;    % If the CT input is too warm, then there is 
end %            no (positive) value of SA that represents frozen seawater. 

if transposed
    SA_freezing = SA_freezing';
end

end