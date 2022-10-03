function travel_time = gsw_travel_time(SA,CT,p,lat)

% gsw_travel_time                vertical acoustic (round-trip) travel time
%                                                        (75-term equation)
%==========================================================================
%
% USAGE:
%  travel_time = gsw_travel_time(SA,CT,p,lat)
%
% DESCRIPTION:
%  Calculates the round-trip acoustic travel time for a path from the
%  bottle concerned up the vertical water column to the sea surface and
%  back to the bottle.
%
%  This function evaluates the pressure integral of specific volume divided
%  by the product of sound speed and the gravitational acceleration, grav,
%  (which is a function of latitude and pressure).
%
%  In order to avoid nonlinear equation of state effects due to the
%  nonlinear dependence of sound speed and specific volume on their input
%  parameters, the vertical data is interpolated so that the data is no
%  more than max_dp_i apart (this is a presure interval). 
% 
%  SA and CT are interpolated using the MRST-PCHIP method of Barker and 
%  McDougall (2020).  This "curve fitting" method uses a Piecewise Cubic 
%  Hermite Interpolating Polynomial to produce a smooth curve with minimal
%  artificial watermasses between the observed data points.  
%
%  The sound speed and specific volume calculations are based on the
%  75-term equations for specific volume (Roquet et al., 2015), as opposed
%  to being based on the Gibbs equations for specific volume and sound 
%  speed).
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  CT    =  in-situ temperature (ITS-90)                          [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  lat  =  latitude in decimal degress north                [ -90 ... +90 ]
%
%  SA & CT need to have the same dimensions.
%  p and lat may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  travel_time  =  vertical acoustic (round-trip) travel time         [ s ]
%
% AUTHOR:
%  Paul Barker, Trevor McDougall and Randy Watts
%                                                      [ help@teos-10.org ]
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
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_travel_time:  Requires four inputs')
end 

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[ml,nl] = size(lat);

if (ms~=mt) | (ns~=nt)
    error('gsw_travel_time: SA & t need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_travel_time: There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    error('gsw_travel_time: need more than one pressure')
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
    error('gsw_travel_time: Inputs array dimensions arguments do not agree')
end 

if (ml == 1) & (nl == 1)              % lat scalar - fill to size of SA
    lat = lat*ones(size(SA));
elseif (ns == nl) & (ml == 1)         % lat is row vector,
    lat = p(ones(1,ms), :);              % copy down each column.
elseif (ms == ml) & (nl == 1)         % lat is column vector,
    lat = lat(:,ones(1,ns));               % copy across each row.
elseif (ns == ml) & (nl == 1)          % lat is a transposed row vector,
    lat = lat.';                              % transposed then
    lat = lat(ones(1,ms), :);                % copy down each column.
elseif (ms == ml) & (ns == nl)
    % ok
else
    error('gsw_travel_time: Inputs array dimensions arguments do not agree')
end 

[Inan] = find(isnan(SA + CT + p));
SA(Inan) = NaN;
CT(Inan) = NaN;
p(Inan) = NaN;

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    lat = lat.';
    transposed = 1;
else
    transposed = 0;
end
[mp,np] = size(p);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%  This max_dp_i is the limit we choose for the evaluation of specific
%  volume in the pressure integration.  That is, the vertical integration
%  of specific volume with respet to pressure is perfomed with the pressure
%  increment being no more than max_dp_i, with the default value being 1
%  dbar.
max_dp_i = 1;
%--------------------------------------------------------------------------

db2Pa = 1e4;

Ishallow = 1:(mp-1);
Ideep = 2:mp;
d_p = (p(Ideep,:) - p(Ishallow,:));

if any(d_p <= 0)
    error('gsw_travel_time: pressure must be monotonic')
end

travel_time = nan(size(SA));
p_surf = 0; % This ensures that there is a bottle at the surface.

%--------------------------------------------------------------------------
% The index [Ibg] (Index-bottle-gaps) indicates where the vertical gaps
% between adjacent "bottles" is greater than max_dp_i.
[Ibg] = find(d_p > max_dp_i);
%--------------------------------------------------------------------------
% The index [Inz] (Index-not-zero) indicates when the shallowest
% "bottle" is not at p = 0 dbar.
[Inz] = find(p(1,:) ~=0);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
if isempty(Ibg) & isempty(Inz)
    % vertical resolution is good (bottle gap is no larger than max_dp_i)
    % & vertical profile begins at the surface (i.e. at p = 0 dbar)
    grav = gsw_grav(lat,p);
        
    v = gsw_specvol(SA,CT,p);
    c = gsw_sound_speed(SA,CT,p);
    B = v./(c.*grav);
    B_av = zeros(size(SA));
    
    %    do the top differently as we do not have a botle pair.
    %top = (2./(sspd(1,:).*dens(1,:).*g(1,:) ) ).* P(1,:)*db2Pascal;
    B_av(1,:) = B(1,:);
    
    B_av(2:mp,:) = 0.5*(B(1:(end-1),:) + B(2:end,:));
    dp = zeros(size(SA));
    dp(2:mp,:) = d_p;
    D = B_av.*dp.*db2Pa;
    
    travel_time = 2.*cumsum(D);
    % The factor of 2 is to allow for the sound to travel both up and down the
    % water column.  "travel_time" is the travel both up and down the water
    % column.  The code will have gotten to here iff the data is "perfect" in
    % the sense that
    %      (i)  it has very fine vertical resolution,
    %     (ii)  each cast starts at p = 0, and
    
else
    % will need to interpolate profiles, doing so one profile at a time.
    for Iprofile = 1:np
        [Inn] = find(~isnan(p(:,Iprofile)));
        % Test if the depth of the cast extends to the surface 
        if (max(p(Inn,Iprofile)) >= p_surf)
            % p_surf is shallower than the pressure of the deepest “bottle” on the
            % vertical profile, thus the dynamic height can be calculated.
            
            % Test if there are vertical gaps between adjacent "bottles" which are
            % greater than max_dp_i, and that there is a "bottle" exactly at the
            % surface pressure.
            [Ibg_i] = find(d_p(:,Iprofile) > max_dp_i);
            if isempty(Ibg_i)
                % Vertical resultion is already good (no larger than max_dp_i, and on this
                % vertical profile there is a "bottle" at exactly p_surf.
                
                %Test if the the shallowest "bottle" is not at p = 0 dbar.
                if min(p(Inn,Iprofile)) > 0
                    %resolution is fine but there is not a bottle at p =0
                    SA_i = SA(Inn(1),Iprofile);
                    SA_i(2:length(Inn)+1) = SA(Inn,Iprofile);
                    CT_i = t(Inn(1),Iprofile);
                    CT_i(2:length(Inn)+1) = CT(Inn,Iprofile);
                    p_i = 0;
                    p_i(2:length(Inn)+1) = p(Inn,Iprofile);
                    [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
                else
                    %resolution is fine and there is a bottle at p =0
                    SA_i = SA(Inn,Iprofile);
                    CT_i = CT(Inn,Iprofile);
                    p_i = p(Inn,Iprofile);
                    [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
                end
            else
                % interpolation is needed.
                p_i = nan(2*round(max(p(Inn,Iprofile)/max_dp_i)),1);
                
                % Test if there is a bottle at p = 0.
                if min(p(Inn,Iprofile)) > 0
                    % There is not a bottle at p = 0.
                    if p_surf < min(p(Inn,Iprofile))
                        dp_iIbottle1 = p_surf;
                        dp_iIbottle2 = p(Inn(1),Iprofile) - p_surf;
                        p_iIbottle1 = 0:dp_iIbottle1/ceil(dp_iIbottle1/max_dp_i):p_surf;
                        p_iIbottle2 = p_surf:dp_iIbottle2/ceil(dp_iIbottle2/max_dp_i):p(Inn(1),Iprofile);
                        p_iIbottle = p_iIbottle1(1:end-1);
                        p_iIbottle((length(p_iIbottle1)):length(p_iIbottle1)+length(p_iIbottle2)-1) =  p_iIbottle2;
                        p_cnt = length(p_iIbottle);
                        p_i(1:p_cnt) = p_iIbottle;
                        top_pad = p_cnt;
                    else
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
                
                % Test for bottle at the surface, if it does not exist then
                % the surface pressure will need to be an interpolated pressure.
                if any(p(Inn,Iprofile) - p_surf == 0)
                    %There is a bottle at the surface. Define interpolation pressures.
                    for Ibottle = 1:(length(Inn)-1)
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        p_iIbottle = p(Inn(Ibottle),Iprofile):dp_iIbottle/ceil(dp_iIbottle/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                        p_cnt_ld = p_cnt+1;
                        p_cnt = p_cnt + length(p_iIbottle(2:length(p_iIbottle)));
                        p_i(p_cnt_ld:p_cnt) = p_iIbottle(2:length(p_iIbottle));
                    end
                else
                    %There is not a bottle at the surface. Define interpolation
                    %pressures to include the surface.
                    for Ibottle = 1:(length(Inn)-1)
                        % Test if the bottle pair spans the reference pressure
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        if (p(Inn(Ibottle+1),Iprofile) - p_surf > 0) & (p(Inn(Ibottle),Iprofile) - p_surf < 0)
                            % The code should not go here as the surface pressure will not be spanned by a bottle pair
                            % surface pressure is spanned by bottle pairs,
                            % need to include the surface as an interpolated
                            % pressure.
                            dp_iIbottle1 = p_surf - p(Inn(Ibottle),Iprofile);
                            dp_iIbottle2 = p(Inn(Ibottle+1),Iprofile) - p_surf;
                            p_iIbottle1 = p(Inn(Ibottle),Iprofile):dp_iIbottle1/ceil(dp_iIbottle1/max_dp_i):p_surf;
                            p_iIbottle2 = p_surf:dp_iIbottle2/ceil(dp_iIbottle2/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                            p_iIbottle = p_iIbottle1(1:end-1);
                            p_iIbottle((length(p_iIbottle1)):length(p_iIbottle1)+length(p_iIbottle2)-1) =  p_iIbottle2;
                        else
                            %surface pressure is not spanned by bottle pairs.
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

%---------------------------------------------------------------------------
% "Cowboy/cowgirl" oceanographers would not action the next 6 lines of
% code.  Instead these "rough & ready" oceanographers would implement the
% one line of code which linearly interpolates.
               [Intrp] = top_pad:length(p_i);
               [SA_i(Intrp), CT_i(Intrp)] =  gsw_SA_CT_interp(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Intrp));
               if any(isnan(SA_i))
                   [Inan] = find(isnan(SA_i));
                   [SA_i(Inan), CT_i(Inan)] = gsw_linear_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan));
               end  
% The linear interpolation below is for use by "cowboy/cowgirl" oceanographers only
% (i.e. those "rough & ready" oceanographers who do not care about accuracy).
%              [SA_i, CT_i] = gsw_linear_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i);
%---------------------------------------------------------------------------
            end
                        
            p_i = p_i(:);
            v = gsw_specvol(SA_i(:),CT_i(:),p_i(:));
            c = gsw_sound_speed(SA_i(:),CT_i(:),p_i(:));
            
            unique_lat = unique(lat(:,Iprofile));
            if length(unique_lat) == 1
                lat_i = ones(length(p_i),1).*unique_lat;% make lat be the same size as p_i
            else
                lat_i = interp1q(p(:,Iprofile),lat(:,Iprofile),p_i(:));
                [Inan] = find(isnan(lat_i));
                lat_i(Inan) = 9.7963;  % default gravitational acceleration
                                       % of 9.7963 m/s^2 (Griffies, 2004)
            end
            grav = gsw_grav(lat_i(:),p_i(:));
            
            B_i = v./(c.*grav);
            B_i_av = 0.5*(B_i(1:(end-1)) + B_i(2:end));
            Da_i = (B_i_av.*diff(p_i).*db2Pa);
            D_i(2:length(p_i(1:end)+1)) = - cumsum(Da_i);
            travel_time(Ibdata,Iprofile) = 2.*D_i(Iidata);
            %     The factor of 2 is to allow for the sound to travel both up and down
            %     the water column.
            
            clear SA_i CT_i p_i
                
        end
    end
end

if transposed
    travel_time = travel_time.';
end 

end
