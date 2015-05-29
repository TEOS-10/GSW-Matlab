function beta_const_CT = gsw_beta_const_CT(SA,t,p)

% gsw_beta_const_CT                          saline contraction coefficient
%                                      at constant Conservative temperature
%==========================================================================
%
% USAGE:  
%  beta_const_CT = gsw_beta_const_CT(SA,t,p)
%
% DESCRIPTION:
%  Calculates the saline (i. e. haline) contraction coefficient of seawater  
%  at constant Conservative Temperature.
%   
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  t   =   in-situ temperature (ITS-90)                           [ deg C ]
%  p   =   sea pressure                                            [ dbar ]
%          (ie. absolute pressure - 10.1325 dbar)
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  beta_const_CT  =  saline contraction coefficient                  [ kg/g ]
%                  at constant Conservative Temperature
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
%    See Eqn. (2.19.3) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_beta_const_CT:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_beta_const_CT: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_beta_const_CT: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    t = t';
    p = p';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

n0 = 0; 
n1 = 1; 
n2 = 2;

pr0 = zeros(size(p)); 

pt0 = gsw_pt0_from_t(SA,t,p);

gp = gsw_gibbs(n0,n0,n1,SA,t,p);

factora = gsw_gibbs(n1,n1,n0,SA,t,p) - gsw_gibbs(n1,n0,n0,SA,pt0,pr0)./(273.15+pt0);

factor = factora./(gp.*gsw_gibbs(n0,n2,n0,SA,t,p));

beta_const_CT = gsw_gibbs(n0,n1,n1,SA,t,p).*factor - ...
                 gsw_gibbs(n1,n0,n1,SA,t,p)./gp;

if transposed
    beta_const_CT = beta_const_CT';
end

end
