function specvol_anom_t_exact = gsw_specvol_anom_standard_t_exact(SA,t,p)

% gsw_specvol_anom_standard_t_exact                 specific volume anomaly
%==========================================================================
%
% USAGE:  
%  specvol_anom_t_exact = gsw_specvol_anom_standard_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, in-situ 
%  temperature and pressure, using the full TEOS-10 Gibbs function.  The 
%  reference value of Absolute Salinity is SSO and the temperature of the
%  reference value is equal to a Conservative Temperature of 0 degrees C. 
% 
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions, 
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
% 
% OUTPUT:
%  specvol_anom_t_exact  =  specific volume anomaly              [ m^3/kg ]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol_anom_standard_t_exact:  Requires three or five inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns )
    error('gsw_specvol_anom_standard_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
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
    error('gsw_specvol_anom_standard_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SSO = gsw_SSO.*ones(size(SA));
pr0 = zeros(size(SA));
CT0 = pr0;
t_zero = gsw_t_from_CT(SSO,CT0,p);

specvol_anom_t_exact = gsw_gibbs(0,0,1,SA,t,p) - gsw_gibbs(0,0,1,SSO,t_zero,p);

if transposed
    specvol_anom_t_exact = specvol_anom_t_exact.';
end

end
