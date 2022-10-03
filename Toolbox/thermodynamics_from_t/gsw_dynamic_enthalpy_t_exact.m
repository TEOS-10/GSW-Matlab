function dynamic_enthalpy_t_exact = gsw_dynamic_enthalpy_t_exact(SA,t,p)

% gsw_dynamic_enthalpy_t_exact                  dyamic enthalpy of seawater
%==========================================================================
%
% USAGE:
%  dynamic_enthalpy_t_exact = gsw_dynamic_enthalpy_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the dynamic enthalpy of seawater from Absolute Salinity, 
%  in situ temperature and pressure.  Dynamic enthalpy was defined by
%  Young (2010) as the difference between enthalpy and potential enthalpy.
%  Note that this function uses the full TEOS-10 Gibbs function (i.e. the 
%  sum of the IAPWS-09 and IAPWS-08 Gibbs functions, see the TEOS-10 
%  Manual, IOC et al. (2010)).   
%    
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  dynamic_enthalpy_t_exact  =  dynamic enthalpy                   [ J/kg ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org 
%
%  Young, W.R., 2010: Dynamic enthalpy, Conservative Temperature, and the
%   seawater Boussinesq approximation. Journal of Physical Oceanography, 
%   40, 394-400.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_dynamic_enthalpy_t_exact: requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(t); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_dynamic_enthalpy_t_exact: SA and CT must have same dimensions')
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
    error('gsw_dynamic_enthalpy_t_exact: Inputs array dimensions arguments do not agree')
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

CT = gsw_CT_from_t(SA,t,p);
dynamic_enthalpy_t_exact = gsw_enthalpy_t_exact(SA,t,p) - gsw_cp0.*CT;

if transposed
    dynamic_enthalpy_t_exact = dynamic_enthalpy_t_exact.';
end

end
