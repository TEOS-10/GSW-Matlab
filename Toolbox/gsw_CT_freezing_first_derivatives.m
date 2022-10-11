function [CTfreezing_SA, CTfreezing_P] = gsw_CT_freezing_first_derivatives(SA,p,saturation_fraction)

% gsw_CT_freezing_first_derivatives       first derivatives of Conservative
%                                     Temperature at which seawater freezes
%==========================================================================
%
% USAGE:
%  [CTfreezing_SA, CTfreezing_P] = ...
%               gsw_CT_freezing_first_derivatives(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the first derivatives of the Conservative Temperature at
%  which seawater freezes, with respect to Absolute Salinity SA and
%  pressure P (in Pa).  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, completely unsaturated) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  CTfreezing_SA = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to Absolute Salinity at
%                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
%
%  CTfreezing_P  = the derivative of the Conservative Temperature at
%                  freezing (ITS-90) with respect to pressure (in Pa) at
%                  fixed Absolute Salinity                         [ K/Pa ]
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_CT_freezing_first_derivatives: Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end

if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_CT_freezing_first_derivatives: saturation fraction MUST be between zero and one.')
end

[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_CT_freezing_first_derivatives: Inputs array dimensions arguments do not agree')
end %if

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
    error('gsw_CT_freezing_first_derivatives: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

tf = gsw_t_freezing(SA,p,saturation_fraction);
[tf_SA, tf_P] = gsw_t_freezing_first_derivatives(SA,p,saturation_fraction);
[CT_SA_wrt_t, CT_T_wrt_t, CT_P_wrt_t] = gsw_CT_first_derivatives_wrt_t_exact(SA,tf,p);

CTfreezing_SA = CT_SA_wrt_t + CT_T_wrt_t.*tf_SA;
CTfreezing_P = CT_P_wrt_t + CT_T_wrt_t.*tf_P;

if transposed
    CTfreezing_SA = CTfreezing_SA.';
    CTfreezing_P = CTfreezing_P.';
end

end