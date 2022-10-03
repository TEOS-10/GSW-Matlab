function cp_ice = gsw_cp_ice(t,p)

% gsw_cp_ice                                  isobaric heat capacity of ice
%==========================================================================
% 
% USAGE:  
%  cp_ice = gsw_cp_ice(t,p)
%
% DESCRIPTION:
%  Calculates the isobaric heat capacity of ice.
%
% INPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where t is MxN.
%
% OUTPUT:
%  cp_ice  =  heat capacity of ice                           [J kg^-1 K^-1]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (1st June, 2021)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_cp_ice:  Requires two inputs')
end

[mt,nt] = size(t);
[mp,np] = size(p);

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(t));
elseif (nt == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,mt), :);              % copy down each column.
elseif (mt == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nt));               % copy across each row.
elseif (nt == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,mt), :);                % copy down each column.
elseif (mt == mp) & (nt == np)
    % ok
else
    error('gsw_cp_ice: Inputs array dimensions arguments do not agree')
end

if mt == 1
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

cp_ice = -(t + gsw_T0).*gsw_gibbs_ice(2,0,t,p);

if transposed
    cp_ice = cp_ice.';
end

end
