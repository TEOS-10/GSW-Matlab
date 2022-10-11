function osmotic_pressure_t_exact = gsw_osmotic_pressure_t_exact(SA,t,pw)

% gsw_osmotic_pressure_t_exact                             osmotic pressure
%==========================================================================
%
% USAGE:
%  osmotic_pressure_t_exact = gsw_osmotic_pressure_t_exact(SA,t,pw)
%
% DESCRIPTION:
%  Calculates the osmotic pressure of seawater.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  pw  =  sea pressure of the pure water side                      [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  pw may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  osmotic_pressure_t_exact  =  osmotic pressure of seawater       [ dbar ]
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
%    See section 3.41 of this TEOS-10 Manual.
%
%  McDougall T.J. and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_osmotic_pressure_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(pw);

if (mt ~= ms | nt ~= ns)
    error('gsw_osmotic_pressure_t_exact: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % pw is a scalar - fill to size of SA
    pw = pw*ones(size(SA));
elseif (ns == np) & (mp == 1)         % pw is row vector,
    pw = pw(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % pw is column vector,
    pw = pw(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % pw is a transposed row vector,
    pw = pw.';                              % transposed then
    pw = pw(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_osmotic_pressure_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    pw = pw.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

db2Pa = 1e4; % conversion factor from dbar to Pa
rec_g2kg = 1e3;% conversion factor from g to kg

gibbs_pure_water = gsw_gibbs(0,0,0,0,t,pw);

p = pw + 235.4684;        % Initial guess of p, in dbar

df_dp = -db2Pa*(gsw_gibbs(0,0,1,SA,t,p) - SA.*gsw_gibbs(1,0,1,SA,t,p)); % Inital guess of df/dp

for Number_of_iterations = 1:2
    p_old = p;
    f = gibbs_pure_water - rec_g2kg.*gsw_chem_potential_water_t_exact(SA,t,p_old);
    p = p_old - f./df_dp; % this is half way through the modified N-R method
    p_mean = 0.5*(p + p_old);
    df_dp = -db2Pa*(gsw_gibbs(0,0,1,SA,t,p_mean) - SA.*gsw_gibbs(1,0,1,SA,t,p_mean)) ;
    p = p_old - f./df_dp;     
end
   
% After two iterations though the modified Newton-Raphson technique 
% (McDougall and Wotherspoon, 2013) the maximum error is 6x10^-12 dbar.

osmotic_pressure_t_exact = p - pw; % osmotic pressure of seawater, in dbar.

if transposed
    osmotic_pressure_t_exact = osmotic_pressure_t_exact.';
end

end
