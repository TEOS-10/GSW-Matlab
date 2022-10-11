function geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,p,p_ref)

% gsw_geo_strf_dyn_height                            dynamic height anomaly
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:  
%  geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,p,p_ref)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume
%  anomaly from the pressure p of the “bottle” to the reference pressure
%  p_ref.
%
%  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
%  to a given reference pressure.  This is the geostrophic streamfunction 
%  for the difference between the horizontal velocity at the pressure 
%  concerned, p, and the horizontal velocity at p_ref.  Dynamic height 
%  anomaly is the geostrophic streamfunction in an isobaric surface.  The 
%  reference values used for the specific volume anomaly are 
%  SSO = 35.16504 g/kg and CT = 0 deg C.  This function calculates 
%  specific volume anomaly using the computationally efficient 
%  expression for specific volume of Roquet et al. (2015). 
%
%  This function evaluates the pressure integral of specific volume using 
%  SA and CT interpolated using the MRST-PCHIP method of Barker and 
%  McDougall (2020).  This "curve fitting" method uses a Piecewise Cubic 
%  Hermite Interpolating Polynomial to produce a smooth curve with minimal
%  artificial watermasses between the observed data points.  
%
%  Note that the 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  For dynamical oceanography we may 
%  take the 75-term rational function expression for specific volume as
%  essentially reflecting the full accuracy of TEOS-10.  The GSW library 
%  function "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to
%  test if some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_ref =  reference pressure                                     [ dbar ]
%           ( i.e. reference absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_ref needs to be a single value, it can have dimensions 1x1 or Mx1 or  
%  1xN or MxN.
%
% OUTPUT:
%  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
%   Note. If p_ref exceeds the pressure of the deepest “bottle” on a 
%     vertical profile, the dynamic height anomaly for each “bottle” 
%     on the whole vertical profile is returned as NaN.
%
% AUTHOR:  
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2020: Two interpolation methods using 
%   multiply-rotated piecewise cubic hermite interpolating polynomials. 
%   J. Atmosph. Ocean. Tech., 37, pp. 605-619. 
%   doi: 10.1175/JTECH-D-19-0211.1. 
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) and section 3.27 of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
% 
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_geo_strf_dyn_height: Requires four inputs')
end 

unique_p_ref = unique(p_ref);
if ~isscalar(unique_p_ref)
    error('gsw_geo_strf_dyn_height: The reference pressure p_ref must be unique')
end
clear p_ref
p_ref = unique_p_ref;

if p_ref < 0
    error('gsw_geo_strf_dyn_height: The reference pressure p_ref must be positive')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_dyn_height: SA & CT need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_geo_strf_dyn_height: There must be at least 2 values')
end

if (mp == 1) & (np == 1)              
    error('gsw_geo_strf_dyn_height: Need more than one pressure')
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
    error('gsw_geo_strf_dyn_height: Inputs array dimensions arguments do not agree')
end

[Inan] = find(isnan(SA + CT + p));
SA(Inan) = NaN;
CT(Inan) = NaN;
p(Inan) = NaN;

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end
[mp,np] = size(p);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

if max(p(:)) < p_ref
    error('gsw_geo_strf_dyn_height: The reference pressure p_ref is deeper than all bottles')
end

%--------------------------------------------------------------------------
%  This max_dp_i is the limit we choose for the evaluation of specific
%  volume in the pressure integration.  That is, the vertical integration
%  of specific volume with respect to pressure is perfomed with the pressure
%  increment being no more than max_dp_i, with the default value being 1
%  dbar.
max_dp_i = 1;
%--------------------------------------------------------------------------

db2Pa = 1e4;

Ishallow = 1:(mp-1);
Ideep = 2:mp;
d_p = (p(Ideep,:) - p(Ishallow,:));

if any(d_p <= 0)
    error('gsw_geo_strf_dyn_height: pressure must be monotonic')
end

geo_strf_dyn_height = nan(size(SA));
 
%--------------------------------------------------------------------------
% The index [Ibg] (Index-bottle-gaps) indicates where the vertical gaps 
% between adjacent "bottles" is greater than max_dp_i.  
[Ibg] = find(d_p > max_dp_i);
%--------------------------------------------------------------------------

% The index [Inz] (Index-not-zero) indicates when the shallowest 
% "bottle" is not at p = 0 dbar.  
[Inz] = find(p(1,:) ~=0);
%--------------------------------------------------------------------------

% The index [Ibp_ref] (Index-bottle-at-p_ref) indicates when there is a 
% cast which does not have a "bottle" at exactly p_ref, and that cast is 
% deeper than p_ref.  
Ibp_ref = [];
Icount = 1;
Iprofile = 1;
Idha = nan(1,np);
test_p_ref = 1;
while test_p_ref == 1
    if Iprofile <= np
        [Ibrp] = find(p(:,Iprofile) == p_ref);
        if isempty(Ibrp) & (p_ref < max(p(:,Iprofile)))
            Ibp_ref = 1;
            test_p_ref = 0;
        elseif ~isempty(Ibrp)
            Idha(Icount) = Iprofile;
            Icount = Icount + 1;
        end
        if Iprofile == np
            if Icount < np
                Idha(Icount:np) = [];
            end
            test_p_ref = 0;
        end
        Iprofile = Iprofile + 1;
    end
end

%--------------------------------------------------------------------------

if isempty(Ibg) & isempty(Inz) & isempty(Ibp_ref)
    % vertical resolution is good (bottle gap is no larger than max_dp_i) 
    % & vertical profile begins at the surface (i.e. at p = 0 dbar) 
    % & each vertical profile contains a "bottle" at exactly p_ref. 
    geo_strf_dyn_height0 = nan(size(SA));
    geo_strf_dyn_height_p_ref = nan(size(SA));
    
    B = gsw_specvol_anom_standard(SA(:,Idha),CT(:,Idha),p(:,Idha));
    
%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and 
% p.  If one wanted to compute dynamic height anomaly with the full TEOS-10
% Gibbs function expression for specific volume, the following lines of 
% code will enable this.
%
%    B = gsw_specvol_anom_standard_CT_exact(SA(:,Idha),CT(:,Idha),p(:,Idha));
%
% Further down the page is a second section which also needs to be
% activated in order to compute dynamic height anomaly with the full 
% TEOS-10 Gibbs function expression for specific volume.
%
%---------------This is the end of the alternative code--------------------

    B_av = zeros(size(SA(:,Idha)));
    B_av(2:mp,:) = 0.5*(B(1:(end-1),:) + B(2:end,:));    
    dp = zeros(size(SA(:,Idha)));
    dp(2:mp,:) = d_p(:,Idha);
    D = B_av.*dp.*db2Pa;
    
    geo_strf_dyn_height0(:,Idha) = -cumsum(D);   
    % "geo_strf_dyn_height0" is the dynamic height anomaly with respect
    % to p_ref = 0 (the surface).  

    geo_strf_dyn_height_p_ref(:,Idha) = meshgrid(geo_strf_dyn_height0(p == p_ref),[1:mp]);
    % "geo_strf_dyn_height_p_ref" is the dynamic height anomaly at p_ref 
    % with respect to the surface.  
   
    geo_strf_dyn_height = geo_strf_dyn_height0 - geo_strf_dyn_height_p_ref;
    % "geo_strf_dyn_height" is the dynamic height anomaly with respect 
    % to p_ref, and is returned.  The code will  have gotten to here iff the 
    % data is "perfect" in the sense that
    %      (i)  it has very fine vertical resolution, 
    %     (ii)  each cast starts at p = 0, and
    %    (iii)  every cast contains a bottle at exactly p_ref.  

else
    % will need to interpolate profiles, doing so one profile at a time.
    for Iprofile = 1:np
        [Inn] = find(~isnan(p(:,Iprofile)));
        % test that there is more than 1 bottle in the profile        
        if isempty(Inn) | length(Inn) == 1
            continue
        end
        % Test if the depth of the cast extends to the reference pressure
        if (max(p(Inn,Iprofile)) >= p_ref)     
        % p_ref is shallower than the pressure of the deepest “bottle” on the
        % vertical profile, thus the dynamic height can be calculated.

        % Test if there are vertical gaps between adjacent "bottles" which are
        % greater than max_dp_i, and that there is a "bottle" exactly at the 
        % reference pressure.
            [Ibg_i] = find(d_p(:,Iprofile) > max_dp_i);
            [Ibrp] = find(p(Inn,Iprofile) == p_ref);
            if isempty(Ibg_i) & ~isempty(Ibrp) 
            % Vertical resultion is already good (no larger than max_dp_i, and on this 
            % vertical profile there is a "bottle" at exactly p_ref. 
                
                % Test if the the shallowest "bottle" is not at p = 0 dbar.  
                if min(p(Inn,Iprofile)) > 0
                    % resolution is fine and there is a bottle at p_ref, but
                    % there is not a bottle at p =0
                    SA_i = SA(Inn(1),Iprofile);
                    SA_i(2:length(Inn)+1) = SA(Inn,Iprofile);
                    CT_i = CT(Inn(1),Iprofile);
                    CT_i(2:length(Inn)+1) = CT(Inn,Iprofile);
                    p_i = 0;
                    p_i(2:length(Inn)+1) = p(Inn,Iprofile);
                    [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
                    [Ibpr] = find(p_i == p_ref);
                else
                    % resolution is fine, there is a bottle at p_ref, and
                    % there is a bottle at p =0
                    SA_i = SA(Inn,Iprofile);
                    CT_i = CT(Inn,Iprofile);
                    p_i = p(Inn,Iprofile);
                    [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
                    [Ibpr] = find(p_i == p_ref);
                end
            else
                % interpolation is needed.
                p_i = nan(2*round(max(p(Inn,Iprofile)/max_dp_i)),1);
               
                % Test if there is a bottle at p = 0.
                if min(p(Inn,Iprofile)) > 0
                    % there is not a bottle at p = 0.
                    % Test if p_ref is shallower than the minimum bottle pressure of the profile
                    if p_ref < min(p(Inn,Iprofile))
                        % p_ref is shallower than the minimum bottle pressure.                       
                        dp_iIbottle1 = p_ref;
                        dp_iIbottle2 = p(Inn(1),Iprofile) - p_ref;
                        p_iIbottle1 = 0:dp_iIbottle1/ceil(dp_iIbottle1/max_dp_i):p_ref;
                        p_iIbottle2 = p_ref:dp_iIbottle2/ceil(dp_iIbottle2/max_dp_i):p(Inn(1),Iprofile);
                        p_iIbottle = p_iIbottle1(1:end-1);
                        p_iIbottle((length(p_iIbottle1)):length(p_iIbottle1)+length(p_iIbottle2)-1) =  p_iIbottle2;
                        p_cnt = length(p_iIbottle);
                        p_i(1:p_cnt) = p_iIbottle;
                        top_pad = p_cnt;  
                    else
                        % p_ref is deeper than the minimum bottle pressure. 
                        p_i(1) = 0;
                        p_i(2) = min(p(Inn,Iprofile));
                        top_pad = 2;
                        p_cnt = 2; 
                    end
                else
                    % there is a bottle at p = 0.
                    p_i(1) = min(p(Inn,Iprofile));
                    top_pad = 1;
                    p_cnt = 1;
                end
                
                % Test for bottle at p_ref, if it does not exist then the reference 
                % pressure will need to be an interpolated pressure.
                if any(p(Inn,Iprofile) - p_ref == 0)
                    % There is a bottle at p_ref. Define interpolation
                    % pressures.
                    for Ibottle = 1:(length(Inn)-1)
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        p_iIbottle = p(Inn(Ibottle),Iprofile):dp_iIbottle/ceil(dp_iIbottle/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                        p_cnt_ld = p_cnt+1;
                        p_cnt = p_cnt + length(p_iIbottle(2:length(p_iIbottle)));
                        p_i(p_cnt_ld:p_cnt) = p_iIbottle(2:length(p_iIbottle));
                    end
                else
                    % There is not a bottle at p_ref. Define interpolation
                    % pressures to include p_ref.
                    for Ibottle = 1:(length(Inn)-1)
                        % Test if the bottle pair spans the reference pressure
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        if (p(Inn(Ibottle+1),Iprofile) - p_ref > 0) & (p(Inn(Ibottle),Iprofile) - p_ref < 0)
                            % reference pressure is spanned by bottle pairs,
                            % need to include p_ref as an interpolated
                            % pressure.
                            dp_iIbottle1 = p_ref - p(Inn(Ibottle),Iprofile);
                            dp_iIbottle2 = p(Inn(Ibottle+1),Iprofile) - p_ref;
                            p_iIbottle1 = p(Inn(Ibottle),Iprofile):dp_iIbottle1/ceil(dp_iIbottle1/max_dp_i):p_ref;
                            p_iIbottle2 = p_ref:dp_iIbottle2/ceil(dp_iIbottle2/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                            p_iIbottle = p_iIbottle1(1:end-1);
                            p_iIbottle((length(p_iIbottle1)):length(p_iIbottle1)+length(p_iIbottle2)-1) =  p_iIbottle2;
                        else
                            % reference pressure is not spanned by bottle pairs.
                            p_iIbottle = p(Inn(Ibottle),Iprofile):dp_iIbottle/ceil(dp_iIbottle/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                        end
                        p_cnt_ld = p_cnt+1;
                        p_cnt = p_cnt + length(p_iIbottle(2:length(p_iIbottle)));
                        p_i(p_cnt_ld:p_cnt) = p_iIbottle(2:length(p_iIbottle));
                    end
                end
                p_i(p_cnt+1:end) = [];
                p_i = p_i(:);
                SA_i = nan(size(p_i));
                CT_i = SA_i;

                [dummy, Iidata, Ibdata] = intersect(p_i,p(:,Iprofile));
                [Ibpr] = find(p_i == p_ref);
                                           
%---------------------------------------------------------------------------
% "Cowboy/cowgirl" oceanographers would not action the next 7 lines of
% code.  Instead these "rough & ready" oceanographers would implement the
% one line of code which linearly interpolates.  
               [Intrp] = top_pad:length(p_i);
               [SA_i(Intrp), CT_i(Intrp)] =  gsw_SA_CT_interp(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Intrp));
               if any(isnan(SA_i))
                   [Inan] = find(isnan(SA_i));
                  % [SA_i(Inan), CT_i(Inan)] = gsw_spline_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan),0.8,10000);
                   [SA_i(Inan), CT_i(Inan)] = gsw_linear_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan));
               end  
                
% The linear interpolation below is for use by "cowboy/cowgirl" oceanographers only 
% (i.e. those "rough & ready" oceanographers who do not care about accuracy).
%              [SA_i, CT_i] = gsw_linear_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i);
%---------------------------------------------------------------------------
            end
            
            p_i = p_i(:);
            B_i = gsw_specvol_anom_standard(SA_i(:),CT_i(:),p_i(:));

%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly using the computationally
% efficient 75-term expression for specific volume in terms of SA, CT and  
% p.  If one wanted to compute dynamic height anomaly with the full TEOS-10 
% Gibbs function expression for specific volume, the following lines of 
% code will enable this.
%
%    B_i = gsw_specvol_anom_standard_CT_exact(SA_i,CT_i,p_i);
%    B_i = B_i(:);
%
%---------------This is the end of the alternative code--------------------
            
            B_i_av = 0.5*(B_i(1:(end-1)) + B_i(2:end));
            D_i = B_i_av.*diff(p_i).*db2Pa;
            D0_i(2:length(B_i_av)+1) = -cumsum(D_i);
           
            geo_strf_dyn_height(Ibdata,Iprofile) = D0_i(Iidata) - D0_i(Ibpr);
            clear SA_i CT_i p_i
        end
    end
end

if transposed
   geo_strf_dyn_height = geo_strf_dyn_height.';
end 

end

