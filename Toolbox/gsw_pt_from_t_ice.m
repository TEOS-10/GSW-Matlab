function pt_ice = gsw_pt_from_t_ice(t,p,p_ref)

% gsw_pt_from_t_ice                            potential temperature of ice
% =========================================================================
%
% USAGE:
%  pt_ice = gsw_pt_from_t_ice(t,p,p_ref)
%
% DESCRIPTION:
%  Calculates potential temperature of ice Ih with the general reference
%  pressure, p_ref, from in-situ temperature, t.
%
%  A faster gsw routine exists if p_ref is indeed zero dbar.  This routine
%  is "gsw_pt0_from_t_ice(t,p)".
%
% INPUT:
%  t  =  in-situ temperature (ITS-90)                             [ deg C ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  p_ref  =  reference pressure                                    [ dbar ]
%  (If reference pressure is not given then it is assumed that reference
%   pressure is zero).
%
%  p & p_ref (if provided) may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where t is MxN.
%
% OUTPUT:
%  pt_ice  =  potential temperature of ice Ih with reference pressure, 
%             p_ref, on the ITS-90 temperature scale              [ deg C ]
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
%    See appendix I of this TEOS-10 Manual.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3)
    error(['gsw_pt_from_t_ice:  Requires either two or three inputs, '...
        'temperature, pressure and (optional) reference pressure'])
end 

if nargin == 2
    % Assume reference pressure is 0 dbar.
    p_ref = zeros(size(t));
end 

[mt,nt] = size(t);
[mp,np] = size(p);
[mpr,npr] = size(p_ref);

if (mp == 1) & (np == 1)              % p scalar - fill to size of t
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
    error('gsw_pt_from_t_ice: Inputs array dimensions arguments do not agree')
end

if (mpr == 1) & (npr == 1)              % p_ref scalar - fill to size of t
    p_ref = p_ref*ones(size(t));
elseif (nt == npr) & (mpr == 1)         % p_ref is row vector,
    p_ref = p_ref(ones(1,mt), :);              % copy down each column.
elseif (mt == mpr) & (npr == 1)         % p_ref is column vector,
    p_ref = p_ref(:,ones(1,nt));               % copy across each row.
elseif (nt == mpr) & (npr == 1)          % p_ref is a transposed row vector,
    p_ref = p_ref.';                              % transposed then
    p_ref = p_ref(ones(1,mt), :);                % copy down each column.
elseif (mt == mpr) & (nt == npr)
    % ok
else
    error('gsw_pt_from_t_ice: Inputs array dimensions arguments do not agree')
end

if mt == 1
    t = t.';
    p = p.';
    p_ref = p_ref.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This is the starting polynomial for pt of ice Ih.
dp = p - p_ref;

p1 = -2.259745637898635e-4;
p2 =  1.486236778150360e-9;
p3 =  6.257869607978536e-12;
p4 = -5.253795281359302e-7;
p5 =  6.752596995671330e-9;
p6 =  2.082992190070936e-11;

pt_ice = t + dp.*(p1 + (p + p_ref).*(p2 + p3*t) ...
         + t.*(p4 + t.*(p5 + p6*t)));

[Isubzero] = find(pt_ice < -gsw_T0);
pt_ice(Isubzero) = -gsw_T0;

pt_ice(pt_ice < -273) = pt_ice(pt_ice < -273) + 0.05; % We have added 0.05
% to the initial estimate of pt_ice at temperatures less than -273 to
% ensure that it is never less than -273.15.

dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice,p_ref);

true_entropy = -gsw_gibbs_ice_part_t(t,p);

for Number_of_iterations = 1:3
    pt_ice_old = pt_ice;
    dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref) - true_entropy;
    pt_ice = pt_ice_old - dentropy./dentropy_dt; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
    ptm_ice = 0.5.*(pt_ice + pt_ice_old);
    dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
    pt_ice = pt_ice_old - dentropy./dentropy_dt; % this is the end of a full iteration of the modified Newton's method
end

if any(pt_ice < -273)
    [Icold] = find(pt_ice < -273);
    
    q1 = -5.849191185294459e-15;
    q2 =  9.330347971181604e-11;
    q3 =  3.415888886921213e-13;
    q4 =  1.064901553161811e-12;
    q5 = -1.454060359158787e-10;
    q6 = -5.323461372791532e-13;
    
    pt_ice(Icold) = t(Icold) + (p(Icold) - p_ref(Icold)).*(q1 + (p(Icold) + p_ref(Icold)).*(q2 + q3.*t(Icold)) ...
        + t(Icold).*(q4 + t(Icold).*(q5 + q6*t(Icold))));        

    dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice(Icold)+0.01,p_ref(Icold));
    % We have added 0.01 to the initial estimate of pt_ice used in the derivative to
    % ensure that it is never less than -273.15 because the derivative approaches zero
    % at absolute zero.
    for Number_of_iterations = 1:3
        pt_ice_old = pt_ice(Icold);
        dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref(Icold)) - true_entropy(Icold);
        pt_ice(Icold) = pt_ice_old - dentropy./dentropy_dt; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
        ptm_ice = 0.5.*(pt_ice(Icold) + pt_ice_old);        
        ptm_ice = ptm_ice + 0.01;    
        % We have added 0.01 to the estimate of ptm_ice for temperatures less than -273 to
        % ensure that they are never less than -273.15 because the derivative approaches zero
        % at absolute zero and the addition of 0.01 degrees C ensures that when we divide
        % by the derivatve in the modified newton routine the function does not blow up.
        dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref(Icold));
        pt_ice(Icold) = pt_ice_old - dentropy./dentropy_dt; % this is the end of a full iteration of the modified Newton's method        
    end
end

if transposed
    pt_ice = pt_ice.';
end

% For temperatures less than -273.1 degrees C the maximum error is less than
% 2x10^-7 degrees C. For temperatures between -273.1 and 273 the maximum error
% is less than 8x10^-8 degrees C, and for temperatures greater than -273 degrees C the
% maximum error is 1.5x10^-12 degrees C.  These errors are over the whole
% ocean depths with both p and pref varying independently between 0 and
% 10,000 dbar, while the in-situ temperature varied independently between
% -273.15 and +2 degrees C.

end
