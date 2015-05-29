function travel_time = gsw_travel_time(SA,t,p,lat)

% gsw_travel_time                vertical acoustic (round-trip) travel time
%==========================================================================
%
% USAGE:  
%  travel_time = gsw_travel_time(SA,t,p,lat)
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
%  more than max_dp_i apart ( this is a presure interval).  This vertical
%  interpoaltion is done on SA and CT (i.e. Absolute Salinity and 
%  Conservative Temperature) since these are conservative variables under
%  the processes of mixing and interpolation.  
%
%  SA and CT “interpolated” with respect to pressure using a scheme based 
%  on the method of Reiniger and Ross (1968).  Our method uses a weighted 
%  mean of (i) values obtained from linear interpolation of the two nearest
%  data points, and (ii) a linear extrapolation of the pairs of data above 
%  and below.  This "curve fitting" method resembles the use of cubic 
%  splines.  
%
%  The sound speed and specific volume calculations are based on the
%  "_t_exact" versions of the GSW software (as opposed to being based on 
%  the 48-term equations for specific volume and sound speed).  
%
% INPUT:
%  SA    =  Absolute Salinity                                      [ g/kg ]
%  t     =  in situ Temperature (ITS-90)                          [ deg C ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  lat  =  latitude in decimal degress north                [ -90 ... +90 ] 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_ref needs to be a single value, it can have dimensions 1x1 or Mx1 or  
%  1xN or MxN.
%
% OUTPUT:
%  travel_time  =  vertical acoustic (round-trip) travel time         [ s ]
%
% AUTHOR:  
%  Paul Barker, Jeff Dunn, Trevor McDougall and Randy Watts 
%                                                      [ help@teos-10.org ]
%
% VERSION NUMBER: 3.04 (10th December, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Reiniger, R. F. and C. K. Ross, 1968: A method of interpolation with
%   application to oceanographic data. Deep-Sea Res. 15, 185-193.
% 
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_travel_time:  Requires four inputs')
end %if

% unique_lat = unique(lat);
% if ~isscalar(unique_lat)
%     error('gsw_travel_time: The latitude lat must be unique')
% end
% clear lat
% lat = unique_lat;


if any(SA < 0)
    error('gsw_travel_time: The Absolute Salinity must be positive!')
end

[ms,ns] = size(SA);
[mt,nt] = size(t);
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
end %if

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
end %if

[Inan] = find(isnan(SA + t + p));
SA(Inan) = NaN;
t(Inan) = NaN;
p(Inan) = NaN;

if ms == 1
    SA = SA.';
    t = t.';
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
        
    specvol = gsw_specvol_t_exact(SA,t,p);
    sound_speed = gsw_sound_speed_t_exact(SA,t,p);
    B = specvol./(sound_speed.*grav);
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
        
% Test if the depth of the cast extends to the reference pressure
        if (max(p(Inn,Iprofile)) >= p_ref)     
% p_ref is shallower than the pressure of the deepest “bottle” on the
% vertical profile, thus the dynamic height can be calculated.

% Test if there are vertical gaps between adjacent "bottles" which are
% greater than max_dp_i, and that there is a "bottle" exactly at the 
% reference pressure.
            [Ibg_i] = find(d_p(:,Iprofile) > max_dp_i);
            [Ibrp] = find(p(Inn,Iprofile) == p_ref);
            if isempty(Ibg_i) & isempty(Ibrp) 
% Vertical resultion is already good (no larger than max_dp_i, and on this 
% vertical profile there is a "bottle" at exactly p_ref. 
                
 %Test if the the shallowest "bottle" is not at p = 0 dbar.  
                if min(p(Inn,Iprofile)) > 0
                    %resolution is fine and there is a bottle at p_ref, but
                    %there is not a bottle at p =0
                    SA_i = SA(Inn(1),Iprofile);
                    SA_i(2:length(Inn)+1) = SA(Inn,Iprofile);
                    CT_i = CT(Inn(1),Iprofile);
                    CT_i(2:length(Inn)+1) = CT(Inn,Iprofile);
                    p_i = 0;
                    p_i(2:length(Inn)+1) = p(Inn,Iprofile);
                    [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
                    [Ibpr] = find(p_i == p_ref);
                else
                    %resolution is fine, there is a bottle at p_ref, and
                    %there is a bottle at p =0
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
                if any(p - p_ref == 0)
                    %There is a bottle at p_ref. Define interpolation
                    %pressures.
                    for Ibottle = 1:(length(Inn)-1)
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        p_iIbottle = p(Inn(Ibottle),Iprofile):dp_iIbottle/ceil(dp_iIbottle/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                        p_cnt_ld = p_cnt+1;
                        p_cnt = p_cnt + length(p_iIbottle(2:length(p_iIbottle)));
                        p_i(p_cnt_ld:p_cnt) = p_iIbottle(2:length(p_iIbottle));
                    end
                else
                    %There is not a bottle at p_ref. Define interpolation
                    %pressures to include p_ref.
                    for Ibottle = 1:(length(Inn)-1)
% Test if the bottle pair spans the reference pressure
                        dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                        if (p(Inn(Ibottle+1),Iprofile) - p_ref > 0) & (p(Inn(Ibottle),Iprofile) - p_ref < 0)
                            %reference pressure is spanned by bottle pairs,
                            %need to include p_ref as an interpolated
                            %pressure.
                            dp_iIbottle1 = p_ref - p(Inn(Ibottle),Iprofile);
                            dp_iIbottle2 = p(Inn(Ibottle+1),Iprofile) - p_ref;
                            p_iIbottle1 = p(Inn(Ibottle),Iprofile):dp_iIbottle1/ceil(dp_iIbottle1/max_dp_i):p_ref;
                            p_iIbottle2 = p_ref:dp_iIbottle2/ceil(dp_iIbottle2/max_dp_i):p(Inn(Ibottle+1),Iprofile);
                            p_iIbottle = p_iIbottle1(1:end-1);
                            p_iIbottle((length(p_iIbottle1)):length(p_iIbottle1)+length(p_iIbottle2)-1) =  p_iIbottle2;
                        else
                            %reference pressure is not spanned by bottle pairs.
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
% "Cowboy/cowgirl" oceanographers would not action the next 7 lines of
% code.  Instead these "rough & ready" oceanographers would implement the
% one line of code which linearly interpolates.  
                [Intrp] = top_pad:length(p_i);
                SA_i(Intrp) = pinterp_from_p(p(:,Iprofile),SA(:,Iprofile),p_i(Intrp));
                CT_i(Intrp) = pinterp_from_p(p(:,Iprofile),CT(:,Iprofile),p_i(Intrp));
                if any(isnan(SA_i))
                    [Inan] = find(isnan(SA_i));
                    [SA_i(Inan), CT_i(Inan)] = gsw_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan));
                end
                
% The linear interpolation below is for use by "cowboy/cowgirl" oceanographers only 
% (i.e. those "rough & ready" oceanographers who do not care about accuracy).
%              [SA_i, CT_i] = gsw_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i);
%---------------------------------------------------------------------------
            end
            
      % make lat be the same size as p_i
            
      p_i = p_i(:); %  What does this line do????            *************
      t_i = gsw_t_from_CT(SA_i(:),CT_i(:),p_i(:));
      SpV = gsw_specvol_t_exact(SA_i(:),t_i(:),p_i(:));
      Sound = gsw_sound_speed_t_exact(SA_i(:),t_i(:),p_i(:));
      
      grav = gsw_grav(SA_i(:),t_i(:),p_i(:));
      B_i = spV./(Sound.*grav);

            
          B_i_av = 0.5*(B_i(1:(end-1)) + B_i(2:end));
          Da_i = (B_i_av.*diff(p_i).*db2Pa);
          D_i(2:length(p_i(2:end)+1) = - cumsum(Da_i);
          geo_travel_time(Ibdata,Iprofile) = 2.*D_i(Iidata);          
%     The factor of 2 is to allow for the sound to travel both up and down 
%     the water column.  
  
         clear SA_i CT_i t_i p_i
            
            
        end
    end
end

if transposed
   travel_time = travel_time.';
end %if

end

%##########################################################################

function [sdat] = pinterp_from_p(odep,obs,sdep)
% pinterp_from_p.
%==========================================================================
% Interpolate values on arbitrary pressures (Designed for bottle data, but 
% fine for 2db CTD data because it handles any gaps safely).
% INPUT:
%  odep  - vector of pressures of the data.
%  obs   - corresponding data values, with nan indicating any gaps.
%  sdep  - pressure to interpolate to.
% OUTPUT:
%  sdat - interpolated values on at the required pressures.
% AUTHOR: Jeff Dunn.  
%==========================================================================

global rr_int_cnt lin_int_cnt dir_sub_cnt r_extrp_cnt;
grad_lim = [];
maxdis = rr_int([],[],sdep);
odep = odep(:);
obs = obs(:);
sdep = sdep(:);
nlvl = length(sdep);
xfn = [0 300 1200 8000];
yfn = [7 15 75 150];
near_lim = interp1(xfn,yfn,sdep);
far_lim = 2*near_lim;
dir_lim = near_lim/5;
sdat = repmat(NaN,nlvl,1);
jj = find(isnan(obs) | isnan(odep) | odep<0 | odep>12000);
if ~isempty(jj)
    obs(jj) = [];
    odep(jj) = [];
end
if ~isempty(obs)
    jj = find((odep(2:end)-odep(1:end-1))<=0);
    if ~isempty(jj)
        obs(jj+1) = [];
        odep(jj+1) = [];
    end
end
ndeps = length(obs);
if ndeps == 0
    % RETURN if no data
    return;
end

if nargin<3 | isempty(maxdis);
    maxdis = 1;
end

if ndeps < 4 | maxdis == -1
    sidx = (1:nlvl)';
else
    % Reiniger & Ross INTERPOLATION (Reiniger and Ross, 1968)
    sdat = rr_int(obs,odep,sdep,1,maxdis);
    sidx = find(isnan(sdat));
    rr_int_cnt = rr_int_cnt + nlvl - length(sidx);
end

% if ~isempty(sidx)  & ndeps >= 2
%     idx = sidx(find(sdep(sidx)>odep(1) & sdep(sidx)<odep(ndeps)));
%     if ~isempty(idx)
%         oidx = interp1(odep,1:ndeps,sdep(idx));
%         dists = [sdep(idx)-odep(floor(oidx)) odep(ceil(oidx))-sdep(idx)];
%         near = min(dists')';
%         far = max(dists')';
%         interp = idx(find(near<near_lim(idx) | far<far_lim(idx)));
%         if ~isempty(interp)
%             sdat(interp) = interp1(odep,obs,sdep(interp));
%             sidx = find(isnan(sdat));
%             lin_int_cnt = lin_int_cnt + length(interp);
%         end
%     end
% end

% if ~isempty(sidx)
%     idx = round(interp1([-99999; odep; 99999],0:ndeps+1,sdep(sidx)));
%     kk = find(abs(odep(idx)-sdep(sidx)) < near_lim(sidx));
%     for jj = kk(:)'
%         sdj = sdep(sidx(jj));
%         odj = odep(idx(jj));
%         x = sdj-odj;
%         new = nan;
%         if abs(x) > 1.5
%             jll = find(abs(odep-sdj) < far_lim(sidx(jj)));
%             if x > 0
%                 jll = flipud(jll);
%             end
%             if length(jll)<2 | max(abs(odep(jll)-odj)) < abs(x)
%                 jll = [];
%             elseif any(abs(diff(odep(jll))) < 1.5)
%                 ll = jll(1);
%                 for mm = jll(2:end)'
%                     if abs(odep(ll(end))-odep(mm))>1.5  & ...
%                             (length(ll) < 4 | abs(odj-odep(mm)) < abs(x))
%                         ll = [ll mm];
%                     end
%                 end
%                 jll = ll;
%             end
%             if length(jll) >= 2
%                 if abs(max(obs(jll))-min(obs(jll)))<.005
%                     new = obs(jll(1));
%                 else
%                     xog = min(odep(jll));
%                     cc = ([ones([length(jll) 1]) odep(jll)-xog]) \ obs(jll);
%                     new = cc(1) + (sdj-xog)*cc(2);
%                 end
%                 r_extrp_cnt = r_extrp_cnt + 1;
%                 if ~isempty(grad_lim)
%                     ofset = abs(obs(idx(jj))-new);
%                     if ofset>abs(x*grad_lim) | ofset>offlim
%                         new = nan;
%                     end
%                 end
%                 sdat(sidx(jj)) = new;
%             end
%         end
%         if isnan(new) & abs(x)<dir_lim(sidx(jj))
%             sdat(sidx(jj)) = obs(idx(jj));
%             dir_sub_cnt = dir_sub_cnt + 1;
%         end
%     end
%     sidx = find(isnan(sdat));
% end

end

%##########################################################################

function [yi,maxdis] = rr_int(y,x,xi,limchk,maxdis)
%==========================================================================
% References:
%  Reiniger RF & Ross CK. 1968.  A method for interpolation with application
%   to oceanographic data.  Deep-Sea Res. 15: 185-193
% AUTHOR: Jeff Dunn  (12th May 1997)  
%==========================================================================

cfrac = 1/15;
coincid_frac = 1/200;
if nargin<4
    limchk = 1;
elseif isempty(limchk)
    limchk = 1;
end
limchk = 0;
xfn = [0 300 1800 8000];
yfn = [50 200 650 1250];
maxdis = interp1(xfn,yfn,xi);
maxdis = [maxdis(:) maxdis(:).*3];
nobs = length(x);
if nobs < 4
    yi = [];
    return;
end
if size(x,1) == 1;  
    x =x';
end
if size(y,1) == 1;  
    y = y';
end
if size(xi,1) == 1;
    xi =xi';
end
tidx = (1:length(xi))';
yi = repmat(NaN,size(tidx));
if x(nobs)>x(1)
    tidx = tidx(find(xi>=x(2) & xi<=x(nobs-1)));
else
    tidx = tidx(find(xi>=x(nobs-1) & xi<=x(2)));
end
if ~isempty(tidx)
    oidx = interp1(x,(1:length(x))',xi(tidx));
    if ~isempty(tidx)
        ridx = round(oidx);
        coincid = find(abs(xi(tidx)-x(ridx)) < (min(maxdis(:))*coincid_frac));
        if ~isempty(coincid)
            yi(tidx(coincid)) = y(ridx(coincid));
            oidx(coincid) = [];
            tidx(coincid) = [];
        end
    end
    rej = find(oidx<2 | oidx>(length(x)-1));
    if ~isempty(rej)
        oidx(rej) = [];
        tidx(rej) = [];
    end
    if ~isempty(tidx)
        fidx = floor(oidx);
        cidx = ceil(oidx);
        inn = [xi(tidx)-x(fidx) x(cidx)-xi(tidx)];
        out = [x(fidx)-x(fidx-1) x(cidx+1)-x(cidx)];
        sumin = inn(:,1) + inn(:,2);
        outtest = abs(out(:,1) + out(:,2) + sumin);
        minsep = abs(min(out')'./sumin);
        rej = find(sumin>maxdis(tidx,1) | outtest>maxdis(tidx,2) | minsep<cfrac);
        if ~isempty(rej)
            tidx(rej) = [];
            inn(rej,:) = [];
            out(rej,:) = [];
            sumin(rej) = [];
            fidx(rej) = [];
            cidx(rej) = [];
        end
    end
    if ~isempty(tidx)
        % Calculate the inner interpolated and 2 outer extrapolated values:
        intp = y(fidx) + ((y(cidx)-y(fidx)).*inn(:,1)./sumin);
        ext1 = y(fidx) + ((y(fidx)-y(fidx-1)).*inn(:,1)./out(:,1));
        ext2 = y(cidx) + ((y(cidx)-y(cidx+1)).*inn(:,2)./out(:,2));
        % Construct the Reiniger&Ross reference curve equation
        % m = the power variable
        m = 1.7;
        top = (abs(ext1-intp).^m).*ext2 + (abs(intp-ext2).^m).*ext1;
        bot = abs(ext1-intp).^m + abs(intp-ext2).^m;
        kk = find(abs(bot)<1E-4);
        if ~isempty(kk)
            yi(tidx(kk)) = intp(kk);
            kk = find(abs(bot)>=1E-4);
            yi(tidx(kk)) = (intp(kk)+(top(kk)./bot(kk)))./2;
        else
            yi(tidx) = (intp+(top./bot))./2;
        end
        gt = y(fidx) > y(cidx);
        yi(tidx) = max([yi(tidx)'; y(fidx+gt)']);
        yi(tidx) = min([yi(tidx)'; y(cidx-gt)']);
    end
end
end

%##########################################################################

