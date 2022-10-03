function wiggliness = gsw_wiggliness(SA,CT,p)

% gsw_wiggliness                                        variation in a cast
%==========================================================================
% 
% USAGE:  
%  wiggliness = gsw_wiggliness(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the how much variation a cast contains by measuring the 
%  change between sucessive bottles on the SA-CT diagram.
%
%                  a.b
%  wiggliness =  -------.(|a|.|b|)^0.5
%                |a|.|b| 
%
%  where a = (SA(i)- SA(i-1)) - (SA(i+1) - SA(i)), 
%  and b = (CT(i)- CT(i-1)) - (CT(i+1) - CT(i)).
%  This calculates the change in angle between consecutive bottle pairs
%  multiplied by the magitude of the change in the property. 
%
% INPUT:  
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions. 
%  p may have dimensions 1x1 or 1xN or MxN, where SA & CT are MxN and N is
%  the number of casts.
%
% OUTPUT:
%  wiggliness  =  wiggliness   (1xN)                              
%
% AUTHOR:  
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_wiggliness:  Requires three inputs')
end 

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_wiggliness: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('gsw_wiggliness: There must be at least 3 bottles')
end

if (mp == 1) & (np == 1)              % p is a scalar - must be two bottles
    error('gsw_wiggliness:  There must be at least 3 bottles')
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_wiggliness: Inputs array dimensions arguments do not agree')
end 

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

[mp,np] = size(p);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

Ishallow = 1:(mp-1);
Ideep = 2:mp;

dp = (p(Ideep,:) - p(Ishallow,:));
if any(dp <= 0)
    error('gsw_wiggliness: pressure must be monotonic')
end

dSA = (SA(Ideep,:) - SA(Ishallow,:));
dCT = (CT(Ideep,:) - CT(Ishallow,:));

a = dSA(1:end-1,:);
b = 4.*dCT(1:end-1,:);
c = dSA(2:end,:);
d = 4.*dCT(2:end,:);

w_part = sqrt(a.^2 + b.^2).*sqrt(c.^2 + d.^2);
w =  (a.*c + b.*d)./w_part;

sqrt_w_part = sqrt(w_part);

w_scaled = abs(acos(w)).*sqrt_w_part;

wiggliness = sum(w_scaled)./(mp-2);

if transposed
    wiggliness = wiggliness.';
end

end