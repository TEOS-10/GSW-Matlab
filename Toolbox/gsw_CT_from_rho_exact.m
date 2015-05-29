function [CT,CT_multiple] = gsw_CT_from_rho_exact(rho,SA,p)

% gsw_t_from_rho_exact                     in-situ temperature from density
% =========================================================================
%
% USAGE:
%  [t,t_multiple] = gsw_t_from_rho_exact(rho,SA,p)
%
% DESCRIPTION:
%  Calculates the in-situ temperature of a seawater sample, for given
%  values of its density, Absolute Salinity and sea pressure (in dbar).
%
% INPUT:
%  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it.
%     That is, it is 'density', not 'density anomaly'.
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  p    =  sea pressure                                            [ dbar ]
%          ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & SA need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & SA are MxN.
%
% OUTPUT:
%  t  =  in-situ temperature                                      [ deg C ]
%  t_multiple  =  in-situ temperature                             [ deg C ]
%    Note that at low salinities, in brackish water, there are two possible
%      temperatures for a single density.  This programme will output both 
%      valid solutions.  To see this second solution the user must call the 
%      programme with two outputs (i.e. [t,t_multiple]), if there is only 
%      one possible solution and the programme has been called with two 
%      outputs the second variable will be set to NaN.
%
% AUTHOR:
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.03 (29th April, 2013)
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

if ~(nargin==3)
    error('gsw_t_from_rho_exact:  Requires three inputs')
end %if

[md,nd] = size(rho);
[ms,ns] = size(SA);
[mp,np] = size(p);

if (ms ~= md | ns ~= nd)
    error('gsw_t_from_rho_exact: rho and SA must have same dimensions')
end

if (mp == 1) & (np == 1)                    % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)               % p is row vector,
    p = p(ones(1,md), :);                   % copy down each column.
elseif (md == mp) & (np == 1)               % p is column vector,
    p = p(:,ones(1,nd));                    % copy across each row.
elseif (nd == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                 % transposed then
    p = p(ones(1,md), :);                   % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_t_from_rho_exact: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    SA = SA.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

[t,t_multiple] = gsw_t_from_rho_exact(rho,SA,p);

CT = gsw_CT_from_t(SA,t,p);
CT_multiple = gsw_CT_from_t(SA,t_multiple,p);

if transposed
    CT = CT.';
    CT_multiple = CT_multiple.';
end

end
