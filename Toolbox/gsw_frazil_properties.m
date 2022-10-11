function [SA_final, CT_final, w_Ih_final] = gsw_frazil_properties(SA_bulk,h_bulk,p)

% gsw_frazil_properties        Absolute Salinity, Conservative Temperature,
%                                       and the ice mass fraction from bulk
%                                       Absolute Salinity and bulk enthalpy
%==========================================================================
%
% USAGE:
%  [SA_final, CT_final, w_Ih_final] = ...
%                                   gsw_frazil_properties(SA_bulk,h_bulk,p)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), w_Ih_final, which results from given values of the bulk
%  Absolute Salinity, SA_bulk, bulk enthalpy, h_bulk, occuring at pressure
%  p.  The final values of Absolute Salinity, SA_final, and Conservative
%  Temperature, CT_final, of the interstitial seawater phase are also
%  returned.  This code assumes that there is no dissolved air in the
%  seawater (that is, saturation_fraction is assumed to be zero
%  throughout the code).
%
%  When the mass fraction w_Ih_final is calculated as being a positive
%  value, the seawater-ice mixture is at thermodynamic equlibrium.  
%
%  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk, 
%  is sufficiently large (i.e. sufficiently "warm") so that there is no ice 
%  present in the final state.  In this case the final state consists of 
%  only seawater rather than being an equlibrium mixture of seawater and 
%  ice which occurs when w_Ih_final is positive.  Note that when 
%  w_Ih_final = 0, the final seawater is not at the freezing temperature. 
%
%  Note that there is another GSW code,
%  gsw_frazil_properties_potential_poly(SA_bulk,h_pot_bulk,p) which
%  treats potential enthalpy as the conservative variable, while, in
%  contrast, the present code treats in situ enthalpy as the conservative
%  variable during the interaction of seawater and ice Ih.
%
% INPUT:
%  SA_bulk =  bulk Absolute Salinity of the seawater and ice mixture
%                                                                  [ g/kg ]
%  h_bulk  =  bulk enthalpy of the seawater and ice mixture        [ J/kg ]
%  p       =  sea pressure                                         [ dbar ]
%             ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA_bulk and h_bulk must have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA_bulk and
%  h_bulk are MxN.
%
% OUTPUT:
%  SA_final    =  Absolute Salinity of the seawater in the final state, 
%                 whether or not any ice is present.               [ g/kg ]
%  CT_final    =  Conservative Temperature of the seawater in the the final
%                 state, whether or not any ice is present.       [ deg C ]
%  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
%                 If this ice mass fraction is positive, the system is at 
%                 thermodynamic equilibrium.  If this ice mass fraction is 
%                 zero there is no ice in the final state which consists 
%                 only of seawater which is warmer than the freezing 
%                 temperature.                                   [unitless]
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
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014:
%   Melting of ice and sea ice into seawater, and frazil ice formation.
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall, T.J., and S.J. Wotherspoon, 2014: A simple modification of
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied
%   Mathematics Letters, 29, 20-25.
%
%  Anonymous, 2014: Modelling the interaction between seawater and frazil
%   ice.  Manuscript, March 2015.  See Eqns. (8) - (15) of this manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_frazil_properties:  Requires three inputs')
end

[ms,ns] = size(SA_bulk);
[mp,np] = size(p);
[mh,nh] = size(h_bulk);

if (mh ~= ms) & (nh ~= ns)
    error('gsw_frazil_properties: The inputs SA_bulk and h_bulk must have the same dimensions')
end

if (mp == 1) & (np == 1)
    p = p*ones(ms,ns);
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.';
    p = p(ones(1,ms), :);
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_frazil_properties: p has wrong dimensions')
end

if ms == 1
    SA_bulk = SA_bulk.';
    h_bulk = h_bulk.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Initial set-up of the calculation
%--------------------------------------------------------------------------
SA_bulk(SA_bulk < 0) = NaN;  % the input value of SA_bulk must not be negative

saturation_fraction = zeros(size(SA_bulk)); % Throughout this code seawater is
                                          % taken to contain no dissolved air.

SA_final = NaN(size(SA_bulk)); 
CT_final = SA_final;
w_Ih_final = SA_final;

%--------------------------------------------------------------------------
% Finding func0
%--------------------------------------------------------------------------
func0 = h_bulk - gsw_enthalpy_CT_exact(SA_bulk,gsw_CT_freezing(SA_bulk,p,saturation_fraction),p);

%--------------------------------------------------------------------------
% When func0 is zero or positive we can immediately calculate the three
% outputs, as the bulk enthalpy, h_bulk, is too large to allow any ice
% at thermodynamic equilibrium.  The result will be (warm) seawater with no
% frazil ice being present.  The three outputs can be set and the rest of
% this code does not need to be performed.
%--------------------------------------------------------------------------
if any(func0 >= 0)
    [J_toowarm] = find(func0 >= 0);
    w_Ih_final(J_toowarm) = 0;
    SA_final(J_toowarm) = SA_bulk(J_toowarm);
    CT_final(J_toowarm) = gsw_CT_from_enthalpy_exact(SA_bulk(J_toowarm),h_bulk(J_toowarm),p(J_toowarm));
end

%--------------------------------------------------------------------------
% Begin to find the solution for those data points that have func0 < 0,
% implying that the output will be a positive ice mass fraction w_Ih_eq.
%--------------------------------------------------------------------------
[I_solve] = find(func0 < 0);
if ~isempty(I_solve)
    SA_bulk_t = SA_bulk(I_solve);
    h_bulk_t = h_bulk(I_solve);
    p_t = p(I_solve);
    saturation_fraction_t = saturation_fraction(I_solve);
    func0_t = func0(I_solve);
    
%--------------------------------------------------------------------------
% Do a quasi-Newton step with a separate polynomial estimate of the
% derivative of func with respect to the ice mass fraction.  This section
% of the code delivers initial values of both w_Ih_t and SA_t to the rest
% of the more formal modified Newtons Method approach of McDougall and
% Wotherspoon (2014).
%--------------------------------------------------------------------------
    num_f = 5.0e-2;
    num_f2 = 6.9e-7;
    num_p = 2.21;
    dfunc_dw_Ih_mean_poly = (3.347814e+05 - num_f.*func0_t.*(1 + num_f2.*func0_t) - num_p.*p_t);
    w_Ih_t = - func0_t./dfunc_dw_Ih_mean_poly;
    w_Ih_t(w_Ih_t > 0.95) = 0.95;
    SA_t = SA_bulk_t./(1 - w_Ih_t);
    SA_t(SA_t < 0 | SA_t > 120) = NaN; 
    
%--------------------------------------------------------------------------
% Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
% fed into the iterative Newton's Method.
%--------------------------------------------------------------------------
    CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
    hf_t = gsw_enthalpy_CT_exact(SA_t,CTf_t,p_t);
    tf_t = gsw_t_freezing(SA_t,p_t,saturation_fraction_t);
    h_Ihf_t = gsw_enthalpy_ice(tf_t,p_t);
    cp_Ih_t = gsw_cp_ice(tf_t,p_t);
    [h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA_t,CTf_t,p_t);
    [CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SA_t,p_t,saturation_fraction_t);
    [tf_SA, dummy] = gsw_t_freezing_first_derivatives(SA_t,p_t,saturation_fraction_t);
    
    dfunc_dw_Ih = hf_t - h_Ihf_t ...
        - SA_t.*(h_hat_SA + h_hat_CT.*CTf_SA + w_Ih_t.*cp_Ih_t.*tf_SA./(1 - w_Ih_t));
        
%--------------------------------------------------------------------------
% Enter the main McDougall-Wotherspoon (2014) modified Newton-Raphson loop
%--------------------------------------------------------------------------
    for number_of_iterations = 1:3
        if number_of_iterations > 1  % On the first iteration CTf_t and h_pot_Ihf_t are both known
            CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
            hf_t = gsw_enthalpy_CT_exact(SA_t,CTf_t,p_t);
        end
        tf_t = gsw_t_freezing(SA_t,p_t,saturation_fraction_t);
        h_Ihf_t = gsw_enthalpy_ice(tf_t,p_t);
        
        func = h_bulk_t - (1 - w_Ih_t).*hf_t - w_Ih_t.*h_Ihf_t; % The function whose zero we seek
        
        w_Ih_old = w_Ih_t;
        w_Ih_t = w_Ih_old - func./dfunc_dw_Ih; % This is half way through the modified Newton's
                                               % method of McDougall and Wotherspoon (2014).
        w_Ih_mean = 0.5.*(w_Ih_t + w_Ih_old);
        w_Ih_mean(w_Ih_mean > 0.9) = NaN; % This ensures that the mass fraction of ice never exceeds 0.9
        
        SA_t = SA_bulk_t./(1 - w_Ih_mean);
        CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
        hf_t = gsw_enthalpy_CT_exact(SA_t,CTf_t,p_t);
        tf_t = gsw_t_freezing(SA_t,p_t,saturation_fraction_t);
        h_Ihf_t = gsw_enthalpy_ice(tf_t,p_t);
        cp_Ih_t = gsw_cp_ice(tf_t,p_t);
        [h_hat_SA, h_hat_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA_t,CTf_t,p_t);
        [CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SA_t,p_t,saturation_fraction_t);
        [tf_SA, dummy] = gsw_t_freezing_first_derivatives(SA_t,p_t,saturation_fraction_t);
        
        dfunc_dw_Ih = hf_t - h_Ihf_t ...
            - SA_t.*(h_hat_SA + h_hat_CT.*CTf_SA + w_Ih_mean.*cp_Ih_t.*tf_SA./(1 - w_Ih_mean));
        
        w_Ih_t = w_Ih_old - func./dfunc_dw_Ih; % This is the end of one full
                                 % iteration of the modified Newton's method
        w_Ih_t(w_Ih_t > 0.9) = NaN;  % This ensures that the mass fraction of ice never exceeds 0.9
        SA_t = SA_bulk_t./(1 - w_Ih_t);
    end
    
    SA_final(I_solve) = SA_t;
    CT_final(I_solve) = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
    w_Ih_final(I_solve) = w_Ih_t;
    
    if any(w_Ih_t < 0)  % This If loop will only trap cases that are smaller        
        [I] = find(w_Ih_t < 0);        % than zero by just machine precision
        SA_final(I_solve(I)) = SA_bulk_t(I);
        CT_final(I_solve(I)) = gsw_CT_from_enthalpy_exact(SA_bulk_t(I),h_bulk(I),p(I));
        w_Ih_final(I_solve(I)) = 0;
    end
    
end

if transposed
    SA_final = SA_final.';
    CT_final = CT_final.';
    w_Ih_final = w_Ih_final.';
end

end
