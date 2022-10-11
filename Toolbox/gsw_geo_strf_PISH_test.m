function geo_strf_PISH_test = gsw_geo_strf_PISH_test(SA,CT,p,p_ref)

db2Pa = 1e4;

[mp,np] = size(p);

Ishallow = 1:(mp-1);
Ideep = 2:mp;
d_p = (p(Ideep,:) - p(Ishallow,:));
 
Idha = 1;

% paul
B_paul = p(:,Idha).*db2Pa.*gsw_specvol_anom_standard(SA(:,Idha),CT(:,Idha),p(:,Idha));

B_av_paul = zeros(size(SA(:,Idha)));
B_av_paul(2:mp,:) = 0.5*(B_paul(1:(end-1),:) + B_paul(2:end,:));
dp = zeros(size(SA(:,Idha)));
dp(2:mp,:) = d_p(:,Idha);
D_paul = B_av_paul.*dp.*db2Pa;

all_PISH_paul(:,Idha) = cumsum(D_paul);


% tom
B_tom = gsw_specvol_anom_standard(SA(:,Idha),CT(:,Idha),p(:,Idha));

B_av_tom = zeros(size(SA(:,Idha)));
B_av_tom(2:mp,:) = 0.5*(B_tom(1:(end-1),:) + B_tom(2:end,:));
dp = zeros(size(SA(:,Idha)));
dp(2:mp,:) = d_p(:,Idha);
D_tom = B_av_tom.*dp.*db2Pa;

geo_strf_dyn_height0_tom(:,Idha) = -cumsum(D_tom);
% "geo_strf_dyn_height0" is the dynamic height anomaly with respect
% to p_ref = 0 (the surface).  

geo_strf_dyn_height_p_ref_tom(:,Idha) = meshgrid(geo_strf_dyn_height0_tom(p == p_ref),[1:mp]);
% "geo_strf_dyn_height_p_ref" is the dynamic height anomaly at p_ref 
% with respect to the surface.  
   
geo_strf_dyn_height_tom = geo_strf_dyn_height0_tom - geo_strf_dyn_height_p_ref_tom;
% "geo_strf_dyn_height" is the dynamic height anomaly with respect 
% to p_ref, and is returned.  

mid_geo_strf_dyn_height_tom = 0.5*(geo_strf_dyn_height_tom(Ishallow) + geo_strf_dyn_height_tom(Ideep));

all_PISH_tom = cumsum(mid_geo_strf_dyn_height_tom.*d_p);

keyboard

end

