function pt0_ice = gsw_pt0_from_t_ice(t,p)

% gsw_pt0_from_t_ice                           potential temperature of ice
%                                       with a reference pressure of 0 dbar
% =========================================================================
%
% USAGE:
%  pt0_ice = gsw_pt0_from_t_ice(t,p)
%
% DESCRIPTION:
%  Calculates potential temperature of ice Ih with a reference pressure of
%  0 dbar, from in-situ temperature, t.
%
% INPUT:
%  t   =  in-situ temperature  (ITS-90)                           [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where t is MxN.
%
% OUTPUT:
%  pt0_ice  =  potential temperature of ice Ih with reference pressure of
%              zero dbar (ITS-90)                                 [ deg C ]
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
%  McDougall T. J. and S. J. Wotherspoon, 2013: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_pt0_from_t_ice:  Requires two inputs')
end

[mt,nt] = size(t);
[mp,np] = size(p);

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
    error('gsw_pt0_from_t_ice: Inputs array dimensions arguments do not agree')
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
    
% This is the starting polynomial for pt0 of ice Ih.
r1 = -2.256611570832386e-4;
r2 = -6.045305921314694e-7;
r3 =  5.546699019612661e-9;
r4 =  1.795030639186685e-11;
r5 =  1.292346094030742e-9;

pt0_ice = t + p.*(r1 + t.*(r2 + t.*(r3 + t.*r4)) + r5*p);
dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

true_entropy = -gsw_gibbs_ice_part_t(t,p);

pt0_ice_old = pt0_ice;
dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
pt0_ice = pt0_ice_old - dentropy./dentropy_dt; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
ptm_ice = 0.5.*(pt0_ice + pt0_ice_old);
dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
pt0_ice = pt0_ice_old - dentropy./dentropy_dt; % this is the end of the first full iteration

if any (t < -45)
    [I_cold] = find(t < -45 & t > -273);
    
    p1 = -2.259745637898635e-4;
    p2 =  1.486236778150360e-9;
    p3 =  6.257869607978536e-12;
    p4 = -5.253795281359302e-7;
    p5 =  6.752596995671330e-9;
    p6 =  2.082992190070936e-11;
    
    pt0_ice(I_cold) = t(I_cold) + p(I_cold).*(p1 + p(I_cold).*(p2 + p3*t(I_cold)) ...
        + t(I_cold).*(p4 + t(I_cold).*(p5 + p6*t(I_cold))));
    
    [Isubzero] = find(pt0_ice(I_cold) < -gsw_T0);
    pt0_ice(I_cold(Isubzero)) = -gsw_T0;
    
    [Isubzero] = find(pt0_ice(I_cold) < -273);
    pt0_ice(I_cold(Isubzero)) = pt0_ice(I_cold(Isubzero)) + 0.05; % We have added 0.05
    % to the initial estimate of pt0_ice at temperatures less than -273 to
    % ensure that it is never less than -273.15.
    
    dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice(I_cold));
        
    for Number_of_iterations = 1:3
        pt0_ice_old = pt0_ice(I_cold);
        dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy(I_cold);
        pt0_ice(I_cold) = pt0_ice_old - dentropy./dentropy_dt; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
        ptm_ice = 0.5.*(pt0_ice(I_cold) + pt0_ice_old);
        dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
        pt0_ice(I_cold) = pt0_ice_old - dentropy./dentropy_dt; % this is the end of a full iteration of the modified Newton's method
    end
    
    if any(pt0_ice < -273)
        [Icold] = find(pt0_ice < -273);
        
        q1 = -5.849191185294459e-15;
        q2 =  9.330347971181604e-11;
        q3 =  3.415888886921213e-13;
        q4 =  1.064901553161811e-12;
        q5 = -1.454060359158787e-10;
        q6 = -5.323461372791532e-13;
        
        pt0_ice(Icold) = t(Icold) + p(Icold).*(q1 + p(Icold).*(q2 + q3.*t(Icold)) ...
            + t(Icold).*(q4 + t(Icold).*(q5 + q6*t(Icold))));
        
        dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice(Icold)+0.01);
        % We have added 0.01 to the initial estimate of pt_ice used in the derivative to
        % ensure that it is never less than -273.15 because the derivative approaches zero
        % at absolute zero.
        for Number_of_iterations = 1:3
            pt0_ice_old = pt0_ice(Icold);
            dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy(Icold);
            pt0_ice(Icold) = pt0_ice_old - dentropy./dentropy_dt; % this is half way through the modified method (McDougall and Wotherspoon, 2013)
            ptm_ice = 0.5.*(pt0_ice(Icold) + pt0_ice_old);
            ptm_ice = ptm_ice + 0.01;
            % We have added 0.01 to the estimate of ptm_ice for temperatures less than -273 to
            % ensure that they are never less than -273.15 because the derivative approaches zero
            % at absolute zero and the addition of 0.01 degrees C ensures that when we divide
            % by the derivatve in the modified newton routine the function does not blow up.
            dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
            pt0_ice(Icold) = pt0_ice_old - dentropy./dentropy_dt; % this is the end of a full iteration of the modified Newton's method
        end
    end
end

if transposed
    pt0_ice = pt0_ice.';
end

% For temperatures less than -273.1 degrees C the maximum error is less than
% 2x10^-7 degrees C. For temperatures between -273.1 and 273 the maximum error
% is less than 8x10^-8 degrees C, and for temperatures greater than -273 degrees C the
% maximum error is 1.5x10^-12 degrees C.   These errors are over the whole
% ocean depths with p varying between 0 and 10,000 dbar, while the in-situ
% temperature varied independently between -273.15 and +2 degrees C.

end
