function [SA_final, CT_final, w_Ih_final] = gsw_frazil_properties_potential(SA_bulk,h_pot_bulk,p)

% gsw_frazil_properties_potential           Absolute Salinity, Conservative
%                               Temperature, and the ice mass fraction when
%                     at thermodynamic equilibrium between seawater and ice 
%==========================================================================
%
% USAGE:
%  [SA_final, CT_final, w_Ih_final] = ...
%                     gsw_frazil_properties_potential(SA_bulk,h_pot_bulk,p)
%
% DESCRIPTION:
%  Calculates the mass fraction of ice (mass of ice divided by mass of ice
%  plus seawater), w_Ih_eq, which results from given values of the bulk
%  Absolute Salinity, SA_bulk, bulk potential enthalpy, h_pot_bulk,
%  occuring at pressure p.  The final equilibrium values of Absolute
%  Salinity, SA_eq, and Conservative Temperature, CT_eq, of the
%  interstitial seawater phase are also returned.  This code assumes that
%  there is no dissolved air in the seawater (that is, saturation_fraction
%  is assumed to be zero thoughout the code).
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
%  Note that this code uses the exact forms of CT_freezing and
%  pot_enthalpy_ice_freezing.
%
% INPUT:
%  SA_bulk     =  bulk Absolute Salinity of the seawater and ice mixture
%                                                                  [ g/kg ]
%  h_pot_bulk  =  bulk potential enthalpy of the seawater and ice mixture
%                                                                  [ J/kg ]
%  p           =  sea pressure                                  [ dbar ]
%                  ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA_bulk and h_pot_bulk must have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA_bulk and
%  h_pot_bulk are MxN.
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
%   ice.  Manuscript, March 2015.  See Eqns. (8)-(15) of this  manuscript.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_frazil_properties_potential:  Requires three inputs')
end

[ms,ns] = size(SA_bulk);
[mp,np] = size(p);
[mh,nh] = size(h_pot_bulk);

if (mh ~= ms) & (nh ~= ns)
    error('gsw_frazil_properties_potential: The inputs SA_bulk and h_pot_bulk must have the same dimensions')
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
    error('gsw_frazil_properties_potential: p has wrong dimensions')
end

if ms == 1
    SA_bulk = SA_bulk.';
    h_pot_bulk = h_pot_bulk.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Initial set-up of the calculation
%--------------------------------------------------------------------------
SA_bulk(SA_bulk < 0) = NaN; % the input value of SA_bulk must not be negative

saturation_fraction = zeros(size(SA_bulk)); % Throughout this code seawater is
                                           % taken to contain no dissolved air
cp0 = gsw_cp0;
rec_cp0 = 2.505092880681252e-04; % rec_cp0 = 1./gsw_cp0

SA_final = NaN(size(SA_bulk));
CT_final = SA_final;
w_Ih_final = SA_final;

%--------------------------------------------------------------------------
% Finding func0.  This is the value of the function, func, that would
% result in the output w_Ih_eq being exactly zero.
%--------------------------------------------------------------------------
func0 = h_pot_bulk - cp0.*gsw_CT_freezing(SA_bulk,p,saturation_fraction);

%--------------------------------------------------------------------------
% Setting the three outputs for data points that have func0 non-negative
%--------------------------------------------------------------------------
if any(func0 >= 0)
    % When func0 is zero or positive then the final answer will contain
    % no frazil ice; that is, it will be pure seawater that is warmer
    % han the freezing temperature.  If func0 >= 0 we do not need to go
    % through the modified Newton-Raphson procedure and we can simply
    % write down the answer, as in the following 4 lines of code.
    [I_toowarm] = find(func0 >= 0);
    w_Ih_final(I_toowarm) = 0;
    SA_final(I_toowarm) = SA_bulk(I_toowarm);
    CT_final(I_toowarm) = rec_cp0.*h_pot_bulk(I_toowarm);
end

%--------------------------------------------------------------------------
% Begin finding the solution for data points that have func0 < 0, so that
% the output will have a positive ice mass fraction w_Ih_eq.
%--------------------------------------------------------------------------
[I_solve] = find(func0 < 0);

if ~isempty(I_solve)
    SA_bulk_t = SA_bulk(I_solve);
    h_pot_bulk_t = h_pot_bulk(I_solve);
    p_t = p(I_solve);
    saturation_fraction_t = saturation_fraction(I_solve);
       
    % Evalaute a polynomial for w_Ih_t in terms of SA_bulk_t, func0 and p_t
    x = SA_bulk_t.*1e-2;  
    y = func0(I_solve).*3.333333333333333e-06;  % y = func0./3e5
    z = p_t.*1e-4;
    
    f01 = -9.041191886754806e-1;
    f02 =  4.169608567309818e-2;
    f03 = -9.325971761333677e-3;
    f04 =  4.699055851002199e-2;
    f05 = -3.086923404061666e-2;
    f06 =  1.057761186019000e-2;
    f07 = -7.349302346007727e-2;
    f08 =  1.444842576424337e-1;
    f09 = -1.408425967872030e-1;
    f10 =  1.070981398567760e-1;
    f11 = -1.768451760854797e-2;
    f12 = -4.013688314067293e-1;
    f13 =  7.209753205388577e-1;
    f14 = -1.807444462285120e-1;
    f15 =  1.362305015808993e-1;
    f16 = -9.500974920072897e-1;
    f17 =  1.192134856624248;
    f18 = -9.191161283559850e-2;
    f19 = -1.008594411490973;
    f20 =  8.020279271484482e-1;
    f21 = -3.930534388853466e-1;
    f22 = -2.026853316399942e-2;
    f23 = -2.722731069001690e-2;
    f24 =  5.032098120548072e-2;
    f25 = -2.354888890484222e-2;
    f26 = -2.454090179215001e-2;
    f27 =  4.125987229048937e-2;
    f28 = -3.533404753585094e-2;
    f29 =  3.766063025852511e-2;
    f30 = -3.358409746243470e-2;
    f31 = -2.242158862056258e-2;
    f32 =  2.102254738058931e-2;
    f33 = -3.048635435546108e-2;
    f34 = -1.996293091714222e-2;
    f35 =  2.577703068234217e-2;
    f36 = -1.292053030649309e-2;
    
    w_Ih_t = y.*(f01 + x.*(f02 + x.*(f03 + x.*(f04 + x.*(f05 + f06*x)))) ...
        + y.*(f07 + x.*(f08 + x.*(f09 + x.*(f10 + f11*x))) + y.*(f12 + x.*(f13 ...
        + x.*(f14 + f15*x)) + y.*(f16 + x.*(f17 + f18*x) + y.*(f19 + f20*x ...
        + f21*y)))) + z.*(f22 + x.*(f23 + x.*(f24 + f25*x)) + y.*(x.*(f26 + f27*x) ...
        + y.*(f28 + f29*x + f30*y)) + z.*(f31 + x.*(f32 + f33*x) + y.*(f34 ...
        + f35.*x + f36*y))));
    
    SA_bulk_t(w_Ih_t > 0.9) = NaN; % The ice mass fraction out of this code
                                        % is restricted to be less than 0.9
    h_pot_bulk_t(isnan(SA_bulk_t)) = NaN;
    w_Ih_t(isnan(SA_bulk_t)) = NaN;
    
    SA_t = SA_bulk_t./(1 - w_Ih_t);          % This is the polynomial-based
                                     % estimate of the Absolute Salinity of
                                          % the interstitial seawater phase
    
%--------------------------------------------------------------------------
% Doing a Newton step with a separate polynomial estimate of the mean
% derivative dfunc_dw_Ih_mean_poly.
%--------------------------------------------------------------------------
    CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
    h_pot_Ihf_t = gsw_pot_enthalpy_ice_freezing(SA_t,p_t);
    func = h_pot_bulk_t - (1 - w_Ih_t).*cp0.*CTf_t - w_Ih_t.*h_pot_Ihf_t;
    
    xa = SA_t.*1e-2;
    
    g01 =  3.332286683867741e5;
    g02 =  1.416532517833479e4;
    g03 = -1.021129089258645e4;
    g04 =  2.356370992641009e4;
    g05 = -8.483432350173174e3;
    g06 =  2.279927781684362e4;
    g07 =  1.506238790315354e4;
    g08 =  4.194030718568807e3;
    g09 = -3.146939594885272e5;
    g10 = -7.549939721380912e4;
    g11 =  2.790535212869292e6;
    g12 =  1.078851928118102e5;
    g13 = -1.062493860205067e7;
    g14 =  2.082909703458225e7;
    g15 = -2.046810820868635e7;
    g16 =  8.039606992745191e6;
    g17 = -2.023984705844567e4;
    g18 =  2.871769638352535e4;
    g19 = -1.444841553038544e4;
    g20 =  2.261532522236573e4;
    g21 = -2.090579366221046e4;
    g22 = -1.128417003723530e4;
    g23 =  3.222965226084112e3;
    g24 = -1.226388046175992e4;
    g25 =  1.506847628109789e4;
    g26 = -4.584670946447444e4;
    g27 =  1.596119496322347e4;
    g28 = -6.338852410446789e4;
    g29 =  8.951570926106525e4;
    
    dfunc_dw_Ih_mean_poly = g01 + xa.*(g02 + xa.*(g03 + xa.*(g04 + g05.*xa))) ...
        + w_Ih_t.*(xa.*(g06 + xa.*(g07 + g08.*xa)) + w_Ih_t.*(xa.*(g09 + g10.*xa) ...
        + w_Ih_t.*xa.*(g11 + g12.*xa + w_Ih_t.*(g13 + w_Ih_t.*(g14 + w_Ih_t.*(g15 ...
        + g16.*w_Ih_t)))))) + z.*(g17 + xa.*(g18 + g19.*xa) + w_Ih_t.*(g20 ...
        + w_Ih_t.*(g21 + g22.*w_Ih_t) + xa.*(g23 + g24.*xa.*w_Ih_t)) ...
        + z.*(g25 + xa.*(g26 + g27.*xa) + w_Ih_t.*(g28 + g29.*w_Ih_t)));
    
    w_Ih_old = w_Ih_t;
    w_Ih_t = w_Ih_old - func./dfunc_dw_Ih_mean_poly; % This is the estimate
                                 % of w_Ih that is fed into Newton's Method
    SA_t = SA_bulk_t./(1 - w_Ih_t);      % This is the initial guess at the 
                          % Absolute Salinity of the interstitial seawater
    
%--------------------------------------------------------------------------
% Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
% fed into Newton's Method.
%--------------------------------------------------------------------------
    CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
    h_pot_Ihf_t = gsw_pot_enthalpy_ice_freezing(SA_t,p_t);
    [CTf_SA, dummy] = gsw_CT_freezing_first_derivatives(SA_t,p_t,saturation_fraction_t);
    [dpot_h_Ihf_dSA, dummy] = gsw_pot_enthalpy_ice_freezing_first_derivatives(SA_t,p_t);
    dfunc_dw_Ih = cp0.*CTf_t - h_pot_Ihf_t - SA_t.*(cp0.*CTf_SA + w_Ih_t.*dpot_h_Ihf_dSA./(1 - w_Ih_t));
    
    if all(w_Ih_t >= 0 & w_Ih_t <= 0.20 & SA_t > 15 & SA_t < 60 & p_t <= 3000)
        max_number_of_iterations = 1;
    elseif all(w_Ih_t >= 0 & w_Ih_t <= 0.85 & SA_t > 0 & SA_t < 120 & p_t <= 3500)
        max_number_of_iterations = 2;
    else
        max_number_of_iterations = 3;
    end
    
    for number_of_iterations = 1:max_number_of_iterations
        if number_of_iterations > 1 % On the first iteration CTf_t and h_pot_Ihf_t are both known
            CTf_t = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
            h_pot_Ihf_t = gsw_pot_enthalpy_ice_freezing(SA_t,p_t);
        end
        func = h_pot_bulk_t - (1 - w_Ih_t).*cp0.*CTf_t - w_Ih_t.*h_pot_Ihf_t; % This is the function, func, whose zero we seek
        w_Ih_old = w_Ih_t;
        w_Ih_t = w_Ih_old - func./dfunc_dw_Ih;
        w_Ih_t(w_Ih_t > 0.9) = NaN;  % This ensures that the mass fraction of ice never exceeds 0.9
        
        SA_t = SA_bulk_t./(1 - w_Ih_t);
    end
    
     SA_final(I_solve) = SA_t;
     CT_final(I_solve) = gsw_CT_freezing(SA_t,p_t,saturation_fraction_t);
     w_Ih_final(I_solve) = w_Ih_t;
    
    if any(w_Ih_t < 0) % This If loop will only trap cases that are smaller 
        [I] = find(w_Ih_t < 0);        % than zero by just machine pecision
        SA_final(I_solve(I)) = SA_bulk_t(I);
        CT_final(I_solve(I)) = rec_cp0.*h_pot_bulk_t(I);
        w_Ih_final(I_solve(I)) = 0;
    end
    
end

if transposed
    SA_final = SA_final.';
    CT_final = CT_final.';
    w_Ih_final = w_Ih_final.';
end

end
