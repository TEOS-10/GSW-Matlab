function SA = gsw_SA_from_rho_t_exact(rho,t,p)

% gsw_SA_from_rho_t_exact       Absolute Salinity from density measurements
%==========================================================================
%
% USAGE:
%  SA = gsw_SA_from_rho_t_exact(rho,t,p)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity of a seawater sample, for given values
%  of its density, in-situ temperature and sea pressure (in dbar). 
%
% INPUT:
%  rho  = density of a seawater sample (e.g. 1026 kg/m^3)        [ kg/m^3 ]
%   Note. This input has not had 1000 kg/m^3 subtracted from it. 
%     That is, it is 'density', not 'density anomaly'.
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  rho & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where rho & t are MxN.
%
% OUTPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%   Note. This is expressed on the Reference-Composition Salinity
%     Scale of Millero et al. (2008). 
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker                      [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.5 of this TEOS-10 Manual. 
%
%  McDougall, T.J., and S.J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin==3)
   error('gsw_SA_from_rho_t_exact:  Requires three inputs')
end %if

[md,nd] = size(rho);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= md | nt ~= nd)
    error('gsw_SA_from_rho_t_exact: rho and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of rho
    p = p*ones(size(rho));
elseif (nd == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,md), :);              % copy down each column.
elseif (md == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,nd));               % copy across each row.
elseif (nd == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,md), :);                          % copy down each column.
elseif (md == mp) & (nd == np)
    % ok
else
    error('gsw_SA_from_rho_t_exact: Inputs array dimensions arguments do not agree')
end %if

if md == 1
    rho = rho.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

v_lab = ones(size(rho))./rho;

v_0 = gsw_gibbs(0,0,1,0,t,p);
v_120 = gsw_gibbs(0,0,1,120,t,p);

SA = 120*(v_lab - v_0)./(v_120 - v_0);            % initial estimate of SA.
SA(SA < 0 | SA > 120) = NaN;

v_SA = (v_120 - v_0)./120; %initial estimate of v_SA, the SA derivative of v

for Number_of_iterations = 1:2
    SA_old = SA;
    delta_v = gsw_gibbs(0,0,1,SA_old,t,p) - v_lab;
    SA = SA_old - delta_v./v_SA ; % this is half way through the modified N-R method (McDougall and Wotherspoon, 2012)
    SA(SA < 0 | SA > 120) = NaN;
    SA_mean = 0.5*(SA + SA_old);
    v_SA = gsw_gibbs(1,0,1,SA_mean,t,p);
    SA = SA_old - delta_v./v_SA;
    SA(SA < 0 | SA > 120) = NaN;
end

% After two iterations of this modified Newton-Raphson iteration,
% the error in SA is no larger than 2x10^-13 g/kg. 

if transposed
    SA = SA.';
end

end
