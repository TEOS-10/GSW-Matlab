function [SP, in_ocean] = gsw_SP_from_SA(SA,p,long,lat)

% gsw_SP_from_SA                  Practical Salinity from Absolute Salinity
%==========================================================================
%
% USAGE:
%  [SP, in_ocean] = gsw_SP_from_SA(SA,p,long,lat)
%
% DESCRIPTION:
%  Calculates Practical Salinity from Absolute Salinity. 
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  p     =  sea pressure                                           [ dbar ]
%          ( ie. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  p, lat and long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SA is MxN.
%
% OUTPUT:
%  SP        =  Practical Salinity  (PSS-78)                   [ unitless ]
%  in_ocean  =  0, if long and lat are a long way from the ocean 
%            =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several 
%    hundred kilometres inland from the coast. 
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd July, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_SP_from_SA:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_SP_from_SA: Inputs array dimensions arguments do not agree')
end %if

[mL,nL] = size(lat);

if (mL == 1) & (nL == 1)             % lat is a scalar - fill to size of SA
    lat = lat*ones(size(SA));
elseif (ns == nL) & (mL == 1)        % lat is a row vector,
    lat = lat(ones(1,ms), :);          % copy down each column.
elseif (ms == mL) & (nL == 1)        % lat is a column vector,
    lat = lat(:,ones(1,ns));           % copy across each row.
elseif (ms == mL) & (ns == nL)
    % ok
else
    error('gsw_SP_from_SA: Inputs array dimensions arguments do not agree')
end %if

[mL,nL] = size(long);
[Iwest] =find(long < 0);
if ~isempty(Iwest)
    long(Iwest) = long(Iwest) + 360; 
end
if (mL == 1) & (nL == 1)            % long is a scalar - fill to size of SA
    long = long*ones(size(SA));
elseif (ns == nL) & (mL == 1)       % long is a row vector,
    long = long(ones(1,ms), :);       % copy down each column.
elseif (ms == mL) & (nL == 1)       % long is a column vector,
    long = long(:,ones(1,ns));        % copy across each row.
elseif (ms == mL) & (ns == nL)
    % ok
else
    error('gsw_SP_from_SA: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    p = p';
    lat = lat';
    long = long';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

inds = find(isfinite(SA)); 

SP = nan(size(SA));
dSA = nan(size(SA));
in_ocean = nan(size(SA));

[dSA(inds), in_ocean(inds)] = gsw_delta_SA(p(inds),long(inds),lat(inds));

SP(inds) = (35/35.16504)*(SA(inds) - dSA(inds));

SP_baltic(inds) = gsw_SP_from_SA_Baltic(SA(inds),long(inds),lat(inds));

indsbaltic = find(~isnan(SP_baltic(inds)));

SP(inds(indsbaltic)) = SP_baltic(inds(indsbaltic));

if transposed
    SP = SP';
    in_ocean = in_ocean';
end

end
