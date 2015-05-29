function result = gsw_delta_SA(p0,longs0,lats0)

%%
% result = gsw_delta_SA(p0,longs0,lats0)
%
% calculate the Absolute Salinity anomaly
%
% p0                  : sea (gauge) pressure               [dbar]
% longs0              : longitude                          [deg E]     
% lats0               : latitude                           [deg N]
%
% result              : Absolute Salinity anomaly          [g/kg]

%%

if gsw_check_arrays(p0,longs0,lats0)
    error('****    input array dimensions in gsw_ASal do not agree    ****')
end

load('gsw_data', 'delta_sa', 'lats', 'longs', 'p', 'ndepth')

nx = length(longs); ny = length(lats); nz = length(p); 

n0 = length(p0);

dlongs = longs(2)-longs(1); dlats = lats(2)-lats(1);

indsx0 = floor(1 + (nx-1)*(longs0-longs(1))./(longs(nx)-longs(1)));
indsx0 = indsx0(:); inds = find(indsx0==nx); indsx0(inds) = nx-1;
              
indsy0 = floor(1 + (ny-1)*(lats0-lats(1))./(lats(ny)-lats(1)));
indsy0 = indsy0(:); inds = find(indsy0==ny); indsy0(inds) = ny-1;

indsz0 = sum(ones(nz,1)*p0(:)' >= p(:)*ones(1,n0));
indsz0 = indsz0(:); 
    
%	adjust in the vertical                                            
                                            
indsn1 = sub2ind([ny,nx], indsy0,     indsx0);      % casts containing data
indsn2 = sub2ind([ny,nx], indsy0,   indsx0+1);
indsn3 = sub2ind([ny,nx], indsy0+1, indsx0+1);
indsn4 = sub2ind([ny,nx], indsy0+1,   indsx0);

nmax = max([ndepth(indsn1)'; ndepth(indsn2)'; ndepth(indsn3)'; ndepth(indsn4)']);

inds1 = find(indsz0(:)' > nmax);                    % casts deeper than GK maximum

p0(inds1) = p(nmax(inds1));

%   have reset p0 here so have to reset indsz0

indsz0 = sum(ones(nz,1)*p0(:)' >= p(:)*ones(1,n0));
indsz0 = indsz0(:); inds = find(indsz0==nz); indsz0(inds) = nz-1;

inds0 = sub2ind([nz,ny,nx], indsz0, indsy0, indsx0);
   
data_indices = [indsx0, indsy0, indsz0, inds0]; data_inds = data_indices(:,3); 
    
r1 = (longs0(:)-longs(indsx0))./(longs(indsx0+1)-longs(indsx0));
s1 = (lats0(:)-lats(indsy0))./(lats(indsy0+1)-lats(indsy0));
t1 = (p0(:)-p(indsz0))./(p(indsz0+1)-p(indsz0));

% rst = [r1,s1,t1]
    
nksum = 0; no_levels_missing = 0;

sa_upper = nan*ones(size(data_inds)); sa_lower = nan*ones(size(data_inds));

sa = nan*ones(size(data_inds));

for k = 1:nz-1
    
  inds_k = find(indsz0==k); nk = length(inds_k);
  
  if nk>0
    nksum = nksum+nk;

    indsx = indsx0(inds_k); indsy = indsy0(inds_k); indsz = k*ones(size(indsx));

    inds_di = find(data_inds==k);

% level k interpolation

    dsa = nan*ones(4,n0);
    inds1 = sub2ind([nz,ny,nx], indsz, indsy,   indsx);
        dsa(1,inds_k) = delta_sa(inds1);      
    inds2 = sub2ind([nz,ny,nx], indsz, indsy,   indsx+1);
        dsa(2,inds_k) = delta_sa(inds2);                             % inds0 + ny*nz
    inds3 = sub2ind([nz,ny,nx], indsz, indsy+1, indsx+1);
        dsa(3,inds_k) = delta_sa(inds3);                             % inds0 + ny*nz + nz
    inds4 = sub2ind([nz ny,nx], indsz, indsy+1,   indsx);
        dsa(4,inds_k) = delta_sa(inds4);                             % inds0 + nz
    
    inds = find(260<=longs0(:)&longs0(:)<=295.217& ...
                    0<=lats0(:)&lats0(:)<=19.55&indsz0(:)==k); 
    if ~isempty(inds)
        dsa(:,inds) = dsa_add_barrier(dsa(:,inds),longs0(inds), ...
            lats0(inds),longs(indsx0(inds)),lats(indsy0(inds)),dlongs,dlats);

    end

    inds = find(isnan(sum(dsa))'&indsz0==k); 
    if ~isempty(inds)
        dsa(:,inds) = dsa_add_mean(dsa(:,inds));
    end
      
    sa_upper(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
                          r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
                            s1(inds_di).*(dsa(4,inds_k)' + ...
                              r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));
          
% level k+1 interpolation

    dsa = nan*ones(4,n0);
    inds1 = sub2ind([nz,ny,nx], indsz+1, indsy,   indsx);
        dsa(1,inds_k) = delta_sa(inds1);      
    inds2 = sub2ind([nz,ny,nx], indsz+1, indsy,   indsx+1);
        dsa(2,inds_k) = delta_sa(inds2);                            % inds1 + ny*nz
    inds3 = sub2ind([nz,ny,nx], indsz+1, indsy+1, indsx+1);
        dsa(3,inds_k) = delta_sa(inds3);                            % inds1 + ny*nz + nz
    inds4 = sub2ind([nz ny,nx], indsz+1, indsy+1,   indsx);
        dsa(4,inds_k) = delta_sa(inds4);                            % inds1 + nz


    inds = find(260<=longs0(:)&longs0(:)<=295.217& ...
                    0<=lats0(:)&lats0(:)<=19.55&indsz0(:)==k); 
    if ~isempty(inds)
        dsa(:,inds) = dsa_add_barrier(dsa(:,inds),longs0(inds), ...
            lats0(inds),longs(indsx0(inds)),lats(indsy0(inds)),dlongs,dlats);
    end
            
    inds = find(isnan(sum(dsa))'&indsz0==k);
    if ~isempty(inds)
        dsa(:,inds) = dsa_add_mean(dsa(:,inds));
    end
            
    sa_lower(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
                          r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
                            s1(inds_di).*(dsa(4,inds_k)' + ...
                              r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));
                        
    inds_different = find(isfinite(sa_upper(inds_di)) & isnan(sa_lower(inds_di)));
            
    if ~isempty(inds_different)
      sa_lower(inds_di(inds_different)) = sa_upper(inds_di(inds_different));
    end

    sa(inds_di) = sa_upper(inds_di) + t1(inds_di).*(sa_lower(inds_di)-sa_upper(inds_di));
            
  else
    no_levels_missing = no_levels_missing+1;
  end
end

inds = find(~isfinite(sa)); sa(inds) = 0;

result = sa;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = dsa_add_mean(dsa)

%%
% result = dsa_add_mean(dsa)
%
% replace nans with namean   
%
% dsa           : absolute salinity anomaly                [g/kg]
%                 of the 4 adjacent neighbours     
%
% result        : nanmean of the 4 adjacent neighbours     [g/kg]

%%

dsa_mean = mean(dsa); inds_nans = find(isnan(dsa_mean)); no_nans = length(inds_nans);

for kk = 1:no_nans
    col = inds_nans(kk);
    inds_kk = find(isnan(dsa(:,col)));
    dsa(inds_kk,col) = nanmean(dsa(:,col));
end

result = dsa;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = dsa_add_barrier(dsa,longs0,lats0,longs,lats,dlongs,dlats)

%% 
% result = dsa_add_barrier(dsa,longs0,lats0,longs,lats,dlongs,dlats)
%
% add a barrier through Central America (Panama) and then average
% over the appropriate side of the barrier
%
% dsa           : absolute salinity anomaly                [g/kg]
% longs0        : longitudes of data                       [deg E]
% lats0         : latitudes of data                        [deg N]
% longs         : longitudes of regular grid               [deg E]
% lats          : latitudes of regular grid                [deg N]
% dlongs         : longitude difference of regular grid    [deg longitude]
% dlats          : latitudes difference of regular grid    [deg latitude]
%
% result        : absolute salinity anomaly of data        [g/kg]

%%

longs_pan = [ 260.0000  272.5900  276.5000  278.6500  280.7300  295.2170];

lats_pan  = [  19.5500   13.9700    9.6000    8.1000    9.3300         0];

lats_lines0 = interp1(longs_pan,lats_pan,longs0);

lats_lines1 = interp1(longs_pan,lats_pan,longs);
lats_lines2 = interp1(longs_pan,lats_pan,longs+dlongs);

for k0 = 1:length(longs0)
    if lats_lines0(k0)<=lats0(k0)
        above_line0 = 1;
    else
        above_line0 = 0;
    end
    if lats_lines1(k0)<=lats(k0)
        above_line(1) = 1;
    else
        above_line(1) = 0;
    end
    if lats_lines1(k0)<=lats(k0)+dlats
        above_line(4) = 1;
    else
        above_line(4) = 0;
    end
    if lats_lines2(k0)<=lats(k0)
        above_line(2) = 1;
    else
        above_line(2) = 0;
    end
    if lats_lines2(k0)<=lats(k0)+dlats
        above_line(3) = 1;
    else
        above_line(3) = 0;
    end
    inds = find(above_line~=above_line0);      % indices of different sides of CA line
    dsa(inds,k0) = nan;
end    
    
dsa_mean = mean(dsa); inds_nans = find(isnan(dsa_mean)); no_nans = length(inds_nans);
            
for kk = 1:no_nans
	col = inds_nans(kk);
	inds_kk = find(isnan(dsa(:,col)));
	dsa(inds_kk,col) = nanmean(dsa(:,col));
end

result = dsa;

end


