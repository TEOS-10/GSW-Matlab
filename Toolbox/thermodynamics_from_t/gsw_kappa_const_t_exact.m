function kappa_const_t_exact = gsw_kappa_const_t_exact(SA,t,p)

% gsw_kappa_const_t_exact                        isothermal compressibility
%==========================================================================
%
% USAGE:  
%  kappa_const_t_exact = gsw_kappa_const_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates isothermal compressibility of seawater. 
%  Note. This is the compressibility of seawater AT CONSTANT IN-SITU
%    TEMPERATURE
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
%  kappa_const_t_exact  =  isothermal compressibility              [ 1/Pa ]
%   Note. The output units are 1/Pa not 1/dbar.
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%   See Eqn. (2.15.1) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_kappa_const_t_exact:  Requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_kappa_const_t_exact: SA and t must have same dimensions')
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
    error('gsw_kappa_const_t_exact: Inputs array dimensions arguments do not agree')
end

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

kappa_const_t_exact = -gsw_gibbs(0,0,2,SA,t,p)./gsw_gibbs(0,0,1,SA,t,p);

if transposed
    kappa_const_t_exact = kappa_const_t_exact.';
end

end