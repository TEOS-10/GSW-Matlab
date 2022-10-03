function dilution_coefficient_t_exact = gsw_dilution_coefficient_t_exact(SA,t,p)

% gsw_dilution_coefficient_t_exact         dilution coefficient of seawater
%==========================================================================
%
% USAGE:
%  dilution_coefficient_t_exact = gsw_dilution_coefficient_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the dilution coefficient of seawater.  The dilution 
%  coefficient of seawater is defined as the Absolute Salinity times the 
%  second derivative of the Gibbs function with respect to Absolute 
%  Salinity, that is, SA.*g_SA_SA.
%
% INPUT:
%  SA =  Absolute Salinity                                         [ g/kg ]
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%        ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  dilution_coefficient_t_exact  =  dilution coefficient   [ (J/kg)(kg/g) ]
%
% AUTHOR: 
%  Trevor McDougall                                    [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
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
   error('gsw_dilution_coefficient_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_dilution_coefficient_t_exact: SA and t must have same dimensions')
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
    error('gsw_dilution_coefficient_t_exact: Inputs array dimensions arguments do not agree')
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

SA(SA < 0) = 0;

sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).

x2 = sfac.*SA;
x = sqrt(x2);
y = t.*0.025;
z = p.*1e-4; %Note.The input pressure (p) is sea pressure in units of dbar.

g08 = 2.0.*(8103.20462414788 + ...
      y.*(2175.341332000392 + y.*(-274.2290036817964 + ...
      y.*(197.4670779425016 + y.*(-68.5590309679152 + 9.98788038278032.*y))) - 90.6734234051316.*z) + ...
      1.5.*x.*(-5458.34205214835 - 980.14153344888.*y + ...
      (4.0./3.0).*x.*(2247.60742726704 - 340.1237483177863.*1.25.*x + 220.542973797483.*y) + ...
      180.142097805543.*z) + ...
      z.*(-219.1676534131548 + (-16.32775915649044 - 120.7020447884644.*z).*z));
    
g08 = x2.*g08 + ... 
      x.*(-7296.43987145382 + z.*(598.378809221703 + ...
          z.*(-156.8822727844005 + (204.1334828179377 - 10.23755797323846.*z).*z)) + ...
          y.*(-1480.222530425046 + z.*(-525.876123559641 + ...
          (249.57717834054571 - 88.449193048287.*z).*z) + ...
          y.*(-129.1994027934126 + z.*(1149.174198007428 + ...
          z.*(-162.5751787551336 + 76.9195462169742.*z)) + ...
          y.*(-30.0682112585625 - 1380.9597954037708.*z + ...
          y.*(2.626801985426835 + 703.695562834065.*z))))) + ...
      11625.62913253464 + 1702.453469893412.*y;
    
dilution_coefficient_t_exact = 0.25.*sfac.*g08;

% Note that this function avoids the singularity that occurs at SA = 0 if
% the straightforward expression for the dilution coefficient of seawater,
% SA*g_SA_SA is simply evaluated as SA.*gsw_gibbs(2,0,0,SA,t,p). 

if transposed
    dilution_coefficient_t_exact = dilution_coefficient_t_exact.';
end

end