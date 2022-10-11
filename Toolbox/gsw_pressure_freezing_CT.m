function pressure_freezing = gsw_pressure_freezing_CT(SA,CT,saturation_fraction)

% gsw_pressure_freezing_CT             pressure of seawater at the freezing
%                                                               temperature
%==========================================================================
%
% USAGE:
%  pressure_freezing = gsw_pressure_freezing_CT(SA,CT,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the pressure (in dbar) of seawater at the freezing
%  temperature.  That is, the output is the pressure at which seawater,
%  with Absolute Salinity SA, Conservative Temperature CT, and with 
%  saturation_fraction of dissolved air, freezes.  If the input values are 
%  such that there is no value of pressure in the range between 0 dbar and 
%  10,000 dbar for which seawater is at the freezing temperature, the 
%  output, pressure_freezing, is put equal to NaN.
%
% INPUT:
%  SA  =  Absolute Salinity of seawater                            [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA and CT are MxN.
%
% OUTPUT:
%  pressure_freezing = sea pressure at which the seawater freezes  [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
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
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of 
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
    error('gsw_pressure_freezing_CT:  Requires either two or three inputs')
end

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
    error('gsw_pressure_freezing_CT: saturation_fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[msf,nsf] = size(saturation_fraction);

if (mt ~= ms | nt ~= ns)
    error('gsw_pressure_freezing_CT: SA and CT must have same dimensions')
end

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_pressure_freezing_CT: Inputs array dimensions arguments do not agree, check saturation_fraction')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0;  % This line ensures that SA is non-negative. 

CT_freezing_p0 = gsw_CT_freezing(SA,zeros(size(SA)),saturation_fraction);
                      %  This is the CT_freezing at zero dbar. 
CT_freezing_p10000 = gsw_CT_freezing(SA,1e4*ones(size(SA)),saturation_fraction);
                      %  This is the CT_freezing at 10,000 dbar.  

% Find CT > CT_freezing_p0.  If this is the case, the input CT value
% represent seawater that will not be frozen at any positive p.  
[Itw] = find(CT > CT_freezing_p0);         % Itw stands for "I_too_warm"
if ~isempty(Itw)
    SA(Itw) = NaN;
    CT(Itw) = NaN;
end
 
% Find CT < CT_freezing_p10000.  If this is the case, the input CT value
% represent seawater that is frozen even at p = 10,000 dbar.   
[Itc] = find(CT < CT_freezing_p10000);         % Itc stands for "I_too_cold"
 if ~isempty(Itc)
    SA(Itc) = NaN; 
    CT(Itc) = NaN;
 end

rec_Pa2dbar = 1e4; 
% The factor rec_Pa2dbar is to have dCTf_dp in units of K/dbar rather than
% K/Pa.  
 
pf = rec_Pa2dbar*(CT_freezing_p0 - CT)./(CT_freezing_p0 - CT_freezing_p10000);
    %  This is the initial (linear) guess of the freezing pressure, in dbar.  

[dummy, CTfreezing_P] = gsw_CT_freezing_first_derivatives(SA,pf,saturation_fraction);
dCTf_dp = rec_Pa2dbar*CTfreezing_P;
    %  This dCTf_dp is the initial value of the partial derivative of 
    %  CT_freezing with respect to pressure (in dbar) at fixed SA, 
    %  assuming that the saturation_fraction is zero. 

Number_of_Iterations = 3;

for I_iter = 1:Number_of_Iterations  
    pf_old = pf;
    f = gsw_CT_freezing(SA,pf_old,saturation_fraction) - CT;
    pf = pf_old - f./dCTf_dp;
                % This is the half-way point of the modified Newton-Raphson
                % method of McDougall and Wotherspoon (2013)
    pfm = 0.5*(pf + pf_old);  % This is now the mean value of p and p_old.
    [dummy, CTfreezing_P] = gsw_CT_freezing_first_derivatives(SA,pfm,saturation_fraction);
    dCTf_dp = rec_Pa2dbar*CTfreezing_P;
    pf = pf_old - f./dCTf_dp; % This is the end of a full iteration of the
                            % modified Newton-Raphson method.   
end

%--------------------------------------------------------------------------
% The following lines of code, if implemented, calculate the error of this 
% function in terms of Conservative Temperature, CT.  
% With Number_of_Iterations = 3, the maximum error in CT is 5x10^-13 C,
% which is equivialnt to a pressure error of 6x10^-10 dbar.  This is the 
% machine precision of the computer.  
%
% SA(SA < 0) = 0;  
% 
% CT_freezing = gsw_CT_freezing(SA,pf,saturation_fraction);
% 
% CT_error = abs(CT_freezing - CT)
% CT_error(p > 10000 | SA > 120 | ...
%      p + SA.*71.428571428571402 > 13571.42857142857) = NaN;
% 
%-----------------This is the end of the error calculation-----------------

pressure_freezing = pf;

%Find any values that are out of the TEOS-10 range, and set them to NaN. 
pressure_freezing(pf > 10000 | SA > 120 | ...
    pf + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if ~isempty(Itw)
    pressure_freezing(Itw) = NaN; % If the CT input is too warm, then there is 
%                no (positive) value of p that represents frozen seawater. 
end

if transposed
    pressure_freezing = pressure_freezing.';
end

end