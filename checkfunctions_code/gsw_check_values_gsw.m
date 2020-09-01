gsw_data = 'gsw_data_chk_cast.mat';
%gsw_data_file = which(['matlab\GSW_extras\',gsw_data]);
gsw_data_file = (gsw_data);

load (gsw_data_file,'SP_chck_cast','t_chck_cast','p_chck_cast',...
    'lat_chck_cast','long_chck_cast');

%load ice test data
load gsw_cv_data.mat
gsw_cv.CT_Arctic = gsw_CT_from_t(gsw_cv.SA_Arctic,gsw_cv.t_Arctic,gsw_cv.p_Arctic);
gsw_cv.t_ice = gsw_cv.t_seaice;
gsw_cv.w_ice = gsw_cv.w_seaice;
gsw_cv.SA_bulk = (1 - gsw_cv.w_ice).*gsw_cv.SA_Arctic;
gsw_cv.h_bulk = (1 - gsw_cv.w_ice).*gsw_enthalpy_CT_exact(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic) + gsw_cv.w_ice.*gsw_enthalpy_ice(gsw_cv.t_ice,gsw_cv.p_Arctic);
gsw_cv.h_pot_bulk = (1 - gsw_cv.w_ice)*gsw_cp0.*gsw_CT_from_t(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic) + gsw_cv.w_ice.*gsw_enthalpy_ice(gsw_cv.t_ice,gsw_cv.p_Arctic);

% gsw_ca_data = 'gsw_chck_vals_errors.mat';
% gsw_data_file = which(gsw_ca_data);
% load (gsw_data_file,'gsw');
SA_chck_cast = gsw_SA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
SA_chck_cast_stable = gsw_stabilise_SA_const_t(SA_chck_cast,t_chck_cast,p_chck_cast,2*gsw_Nsquared_JMcD_lowerlimit(p_chck_cast,lat_chck_cast));
SP_chck_cast = gsw_SP_from_SA(SA_chck_cast_stable,p_chck_cast,long_chck_cast,lat_chck_cast);
clear SA_chck_cast SA_chck_cast_stable
% SA_chck_cast = gsw_SA_from_SP(gsw_cv.SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);

gsw_cv.SP_chck_cast = SP_chck_cast;
gsw_cv.t_chck_cast = t_chck_cast;
gsw_cv.p_chck_cast = p_chck_cast;
gsw_cv.lat_chck_cast = lat_chck_cast;
gsw_cv.long_chck_cast = long_chck_cast;

pr = 0;
gsw_cv.pr = pr;

pr_05 = 5000;
gsw_cv.pr_05 = pr_05;

p_chck_cast_shallow = p_chck_cast;
gsw_cv.p_chck_cast_shallow = p_chck_cast_shallow;

p_chck_cast_deep = p_chck_cast + 10;
gsw_cv.p_chck_cast_deep = p_chck_cast_deep;

delta_p_chck_cast = nan(size(SP_chck_cast));
delta_p_chck_cast(1,:) = p_chck_cast(1,:) - 0;
for I = 2:45
    delta_p_chck_cast(I,:) = p_chck_cast(I,:) - p_chck_cast(I-1,:);
end
gsw_cv.delta_p_chck_cast = delta_p_chck_cast;

Neutral_Density = 26.8*ones(1,length(lat_chck_cast));
gsw_cv.Neutral_Density = Neutral_Density;

p_Neutral_Density = [200 550 91];
gsw_cv.p_Neutral_Density = p_Neutral_Density;

SA_ref = 30;
gsw_cv.SA_ref = SA_ref;

CT_ref = 5;
gsw_cv.CT_ref = CT_ref;

p_i = ([0:100:5000]).';
gsw_cv.p_i = p_i;

%% Practical Salinity (SP):- PSS-78 

C = gsw_C_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.C_from_SP = C;

gsw_cv.SP_from_C = gsw_SP_from_C(C,t_chck_cast,p_chck_cast);

R = gsw_R_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.R_from_SP = R;

gsw_cv.SP_from_R = gsw_SP_from_R(R,t_chck_cast,p_chck_cast);

gsw_cv.Rt_chck_cast = gsw_Rt_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.SP_salinometer = gsw_SP_salinometer(gsw_cv.Rt_chck_cast,t_chck_cast);

gsw_cv.SK_chck_cast = SP_chck_cast;
gsw_cv.SP_from_SK = gsw_SP_from_SK(gsw_cv.SK_chck_cast);

%% Absolute Salinity (SA) and Preformed Salinity (Sstar) 

gsw_cv.SA_chck_cast = gsw_SA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
SA_chck_cast = gsw_cv.SA_chck_cast;
%gsw_cv.SA_from_SP = gsw_cv.SA_chck_cast;
 
gsw_cv.SA_from_SP = gsw_SA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
 
Sstar = gsw_Sstar_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
gsw_cv.Sstar_from_SP = Sstar;
  
%% Conservative Temperature (CT) 

gsw_cv.CT_chck_cast = gsw_CT_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.CT_from_t = gsw_cv.CT_chck_cast ;
CT_chck_cast = gsw_cv.CT_chck_cast ;

%% other conversions between temperatures, salinities, pressure and height 

gsw_cv.deltaSA_from_SP = gsw_deltaSA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);

[gsw_cv.SA_SA_Sstar_from_SP, gsw_cv.Sstar_SA_Sstar_from_SP] = gsw_SA_Sstar_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);

SR = gsw_SR_from_SP(gsw_cv.SP_chck_cast);
gsw_cv.SR_from_SP = SR;

gsw_cv.SP_from_SR = gsw_SP_from_SR(SR);
 
SP_from_SA = gsw_SP_from_SA(SA_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
gsw_cv.SP_from_SA = SP_from_SA;

gsw_cv.Sstar_from_SA = gsw_Sstar_from_SA(SA_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
 
gsw_cv.SA_from_Sstar = gsw_SA_from_Sstar(Sstar,p_chck_cast,long_chck_cast,lat_chck_cast);
 
gsw_cv.SP_from_Sstar = gsw_SP_from_Sstar(Sstar,p_chck_cast,long_chck_cast,lat_chck_cast);

gsw_cv.t_from_CT = gsw_t_from_CT(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
pt = gsw_pt_from_t(SA_chck_cast,t_chck_cast,p_chck_cast,pr);
gsw_cv.pt_from_t = pt;
 
pt0 = gsw_pt0_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.pt0_from_t = pt0;
 
gsw_cv.pt_from_CT = gsw_pt_from_CT(SA_chck_cast,CT_chck_cast);
 
gsw_cv.CT_from_pt = gsw_CT_from_pt(SA_chck_cast,pt);
 
gsw_cv.pot_enthalpy_from_pt = gsw_pot_enthalpy_from_pt(SA_chck_cast,pt);
 
gsw_cv.t90_from_t68 = gsw_t90_from_t68(t_chck_cast);
 
gsw_cv.t90_from_t48 = gsw_t90_from_t48(t_chck_cast);

z_from_p = gsw_z_from_p(p_chck_cast,lat_chck_cast);
gsw_cv.z_from_p = z_from_p;

gsw_cv.p_from_z = gsw_p_from_z(z_from_p,lat_chck_cast);
 
depth_from_z = gsw_depth_from_z(z_from_p);
gsw_cv.depth_from_z = depth_from_z;

z_from_depth = gsw_z_from_depth(depth_from_z);
gsw_cv.z_from_depth = z_from_depth;

Abs_Pressure = gsw_Abs_Pressure_from_p(p_chck_cast);
gsw_cv.Abs_Pressure_from_p = Abs_Pressure;

gsw_cv.p_from_Abs_Pressure = gsw_p_from_Abs_Pressure(Abs_Pressure);

entropy_from_CT = gsw_entropy_from_CT(SA_chck_cast,CT_chck_cast);
gsw_cv.entropy_from_CT = entropy_from_CT;

gsw_cv.CT_from_entropy = gsw_CT_from_entropy(SA_chck_cast,entropy_from_CT);

entropy_from_pt = gsw_entropy_from_pt(SA_chck_cast,pt);
gsw_cv.entropy_from_pt = entropy_from_pt;

gsw_cv.pt_from_entropy = gsw_pt_from_entropy(SA_chck_cast,entropy_from_pt);

entropy_from_t = gsw_entropy_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.entropy_from_t = entropy_from_t;

gsw_cv.t_from_entropy = gsw_t_from_entropy(SA_chck_cast,entropy_from_t,p_chck_cast);

gsw_cv.adiabatic_lapse_rate_from_CT = gsw_adiabatic_lapse_rate_from_CT(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.adiabatic_lapse_rate_from_t = gsw_adiabatic_lapse_rate_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.molality_from_SA = gsw_molality_from_SA(SA_chck_cast);
 
gsw_cv.ionic_strength_from_SA = gsw_ionic_strength_from_SA(SA_chck_cast);

%% specific volume, density and enthalpy 

gsw_cv.specvol = gsw_specvol(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.alpha = gsw_alpha(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.beta = gsw_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.alpha_on_beta = gsw_alpha_on_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_vab, gsw_cv.alpha_vab, gsw_cv.beta_vab] = gsw_specvol_alpha_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
[gsw_cv.v_SA, gsw_cv.v_CT, gsw_cv.v_P] = gsw_specvol_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_SA, gsw_cv.v_SA_CT, gsw_cv.v_CT_CT, gsw_cv.v_SA_P, gsw_cv.v_CT_P] = gsw_specvol_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_wrt_h, gsw_cv.v_h] = gsw_specvol_first_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_SA_wrt_h, gsw_cv.v_SA_h, gsw_cv.v_h_h] = gsw_specvol_second_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.specvol_anom = gsw_specvol_anom(SA_chck_cast,CT_chck_cast,p_chck_cast,SA_ref,CT_ref);
 
gsw_cv.specvol_anom_standard = gsw_specvol_anom_standard(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
rho = gsw_rho(SA_chck_cast,CT_chck_cast,p_chck_cast);
 gsw_cv.rho = rho;

[gsw_cv.rho_rab, gsw_cv.alpha_rab, gsw_cv.beta_rab] = gsw_rho_alpha_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA, gsw_cv.rho_CT, gsw_cv.rho_P] = gsw_rho_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_SA, gsw_cv.rho_SA_CT, gsw_cv.rho_CT_CT, gsw_cv.rho_SA_P, gsw_cv.rho_CT_P] = gsw_rho_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_wrt_h, gsw_cv.rho_h] = gsw_rho_first_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_SA_wrt_h, gsw_cv.rho_SA_h, gsw_cv.rho_h_h] = gsw_rho_second_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.sigma0 = gsw_sigma0(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma1 = gsw_sigma1(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma2 = gsw_sigma2(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma3 = gsw_sigma3(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma4 = gsw_sigma4(SA_chck_cast,CT_chck_cast);

gsw_cv.cabbeling = gsw_cabbeling(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
gsw_cv.thermobaric = gsw_thermobaric(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.enthalpy = gsw_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.enthalpy_diff =  gsw_enthalpy_diff(SA_chck_cast,CT_chck_cast,p_chck_cast_shallow,p_chck_cast_deep);
 
gsw_cv.dynamic_enthalpy = gsw_dynamic_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
[gsw_cv.h_SA, gsw_cv.h_CT] = gsw_enthalpy_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);  
 
[gsw_cv.h_SA_SA, gsw_cv.h_SA_CT, gsw_cv.h_CT_CT] = gsw_enthalpy_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast); 

gsw_cv.sound_speed = gsw_sound_speed(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.kappa = gsw_kappa(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.internal_energy = gsw_internal_energy(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.u_SA, gsw_cv.u_CT, gsw_cv.u_P] = gsw_internal_energy_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.u_SA_SA, gsw_cv.u_SA_CT, gsw_cv.u_CT_CT, gsw_cv.u_SA_P, gsw_cv.u_CT_P] = gsw_internal_energy_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.CT_from_enthalpy = gsw_CT_from_enthalpy(SA_chck_cast,gsw_cv.enthalpy,p_chck_cast);

gsw_cv.SA_from_rho = gsw_SA_from_rho(rho,CT_chck_cast,p_chck_cast);

gsw_cv.CT_from_rho = gsw_CT_from_rho(rho,SA_chck_cast,p_chck_cast);

gsw_cv.CT_maxdensity = gsw_CT_maxdensity(SA_chck_cast,p_chck_cast);

%% vertical stability and interpolation
 
[gsw_cv.Tu, gsw_cv.Rsubrho, gsw_cv.p_mid_TuRsr] = gsw_Turner_Rsubrho(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.n2, gsw_cv.p_mid_n2] = gsw_Nsquared(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast);

[gsw_cv.n2min, gsw_cv.n2min_pmid, gsw_cv.n2min_specvol, gsw_cv.n2min_alpha, gsw_cv.n2min_beta, gsw_cv.n2min_dsa, gsw_cv.n2min_dct, gsw_cv.n2min_dp] = gsw_Nsquared_min(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast);

gsw_cv.mlp = gsw_mlp(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.n2_lowerlimit = gsw_Nsquared_lowerlimit(p_chck_cast,long_chck_cast,lat_chck_cast);

[gsw_cv.SAi_SACTinterp, gsw_cv.CTi_SACTinterp] = gsw_SA_CT_interp(SA_chck_cast,CT_chck_cast,p_chck_cast,p_i);

gsw_cv.ti_tinterp = gsw_t_interp(t_chck_cast,p_chck_cast,p_i);

[gsw_cv.traceri_tracerCTinterp, gsw_cv.CTi_tracerCTinterp] = gsw_tracer_CT_interp(SA_chck_cast,CT_chck_cast,p_chck_cast,p_i,9);

gsw_cv.traceri_tracerinterp = gsw_tracer_interp(CT_chck_cast,p_chck_cast,p_i);

[gsw_cv.IPVfN2, gsw_cv.p_mid_IPVfN2] = gsw_IPV_vs_fNsquared_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
  
%% geostrophic streamfunctions, travel time and Geostrophic velocity %%- Tomlab optimisation

geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
gsw_cv.geo_strf_dyn_height = geo_strf_dyn_height;
 
[gsw_cv.geo_strf_dyn_height_pc, gsw_cv.geo_strf_dyn_height_pc_p_mid] = gsw_geo_strf_dyn_height_pc(SA_chck_cast,CT_chck_cast,delta_p_chck_cast);
 
gsw_cv.geo_strf_isopycnal = gsw_geo_strf_isopycnal(SA_chck_cast,CT_chck_cast,p_chck_cast,pr,Neutral_Density,p_Neutral_Density);
 
[gsw_cv.geo_strf_isopycnal_pc, gsw_cv.geo_strf_isopycnal_pc_p_mid] = gsw_geo_strf_isopycnal_pc(SA_chck_cast,CT_chck_cast,delta_p_chck_cast,26.8,3);
  
gsw_cv.geo_strf_Montgomery = gsw_geo_strf_Montgomery(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
 
gsw_cv.geo_strf_Cunningham = gsw_geo_strf_Cunningham(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);

gsw_cv.geo_strf_steric_height = gsw_geo_strf_steric_height(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);

gsw_cv.geo_strf_PISH = gsw_geo_strf_PISH(SA_chck_cast,CT_chck_cast,p_chck_cast,pr_05);

gsw_cv.travel_time = gsw_travel_time(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast(1));

[gsw_cv.geo_strf_velocity, gsw_cv.geo_strf_velocity_mid_lat, gsw_cv.geo_strf_velocity_mid_long] = gsw_geostrophic_velocity(geo_strf_dyn_height,long_chck_cast,lat_chck_cast,p_chck_cast);

%% geostrophic streamfunctions, travel time and Geostrophic velocity - IBM optimisation

% geo_strf_dyn_height_IBMoptim = gsw_geo_strf_dyn_height_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% gsw_cv.geo_strf_dyn_height_IBMoptim = geo_strf_dyn_height_IBMoptim;
%   
% gsw_cv.geo_strf_isopycnal_IBMoptim = gsw_geo_strf_isopycnal_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr,Neutral_Density,p_Neutral_Density);
%    
% gsw_cv.geo_strf_Montgomery_IBMoptim = gsw_geo_strf_Montgomery_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
%  
% gsw_cv.geo_strf_Cunningham_IBMoptim = gsw_geo_strf_Cunningham_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_steric_height_IBMoptim = gsw_geo_strf_steric_height_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_PISH_IBMoptim = gsw_geo_strf_PISH_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr_05);
% 
% gsw_cv.travel_time_IBMoptim = gsw_travel_time_IBMoptim(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast(1));
% 
% %% geostrophic streamfunctions, travel time and Geostrophic velocity - matlab optimisation
% 
% geo_strf_dyn_height_optim = gsw_geo_strf_dyn_height_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% gsw_cv.geo_strf_dyn_height_optim = geo_strf_dyn_height_optim;
%   
% gsw_cv.geo_strf_isopycnal_optim = gsw_geo_strf_isopycnal_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr,Neutral_Density,p_Neutral_Density);
%    
% gsw_cv.geo_strf_Montgomery_optim = gsw_geo_strf_Montgomery_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
%  
% gsw_cv.geo_strf_Cunningham_optim = gsw_geo_strf_Cunningham_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_steric_height_optim = gsw_geo_strf_steric_height_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_PISH_optim = gsw_geo_strf_PISH_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr_05);
% 
% gsw_cv.travel_time_optim = gsw_travel_time_optim(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast(1));
% 
% %% geostrophic streamfunctions, travel time and Geostrophic velocity - no optimisation
% 
% geo_strf_dyn_height_nooptim = gsw_geo_strf_dyn_height_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% gsw_cv.geo_strf_dyn_height_nooptim = geo_strf_dyn_height_nooptim;
%   
% gsw_cv.geo_strf_isopycnal_nooptim = gsw_geo_strf_isopycnal_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr,Neutral_Density,p_Neutral_Density);
%    
% gsw_cv.geo_strf_Montgomery_nooptim = gsw_geo_strf_Montgomery_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
%  
% gsw_cv.geo_strf_Cunningham_nooptim = gsw_geo_strf_Cunningham_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_steric_height_nooptim = gsw_geo_strf_steric_height_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
% 
% gsw_cv.geo_strf_PISH_nooptim = gsw_geo_strf_PISH_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,pr_05);
% 
% gsw_cv.travel_time_nooptim = gsw_travel_time_nooptim(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast(1));

%% neutral versus isopycnal slopes and ratios 

gsw_cv.isopycnal_slope_ratio = gsw_isopycnal_slope_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
 
[gsw_cv.G_CT, gsw_cv.p_mid_G_CT] = gsw_isopycnal_vs_ntp_CT_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
 
gsw_cv.ntpptCT = gsw_ntp_pt_vs_CT_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast);

%% derivatives of entropy, CT and pt

[gsw_cv.CT_SA, gsw_cv.CT_pt] = gsw_CT_first_derivatives(SA_chck_cast,pt);  
 
[gsw_cv.CT_SA_SA, gsw_cv.CT_SA_pt, gsw_cv.CT_pt_pt] = gsw_CT_second_derivatives(SA_chck_cast,pt);
  
[gsw_cv.eta_SA, gsw_cv.eta_CT] = gsw_entropy_first_derivatives(SA_chck_cast,CT_chck_cast);
 
[gsw_cv.eta_SA_SA, gsw_cv.eta_SA_CT, gsw_cv.eta_CT_CT] = gsw_entropy_second_derivatives(SA_chck_cast,CT_chck_cast);
 
[gsw_cv.pt_SA, gsw_cv.pt_CT] = gsw_pt_first_derivatives(SA_chck_cast,CT_chck_cast);
 
[gsw_cv.pt_SA_SA, gsw_cv.pt_SA_CT, gsw_cv.pt_CT_CT] = gsw_pt_second_derivatives(SA_chck_cast,CT_chck_cast);  

%% seawater properties at freezing temperatures

gsw_cv.CT_freezing = gsw_CT_freezing(SA_chck_cast,p_chck_cast,0.5);

gsw_cv.CT_freezing_poly = gsw_CT_freezing_poly(SA_chck_cast,p_chck_cast,0.5);

gsw_cv.t_freezing = gsw_t_freezing(SA_chck_cast,p_chck_cast,0.5);

gsw_cv.t_freezing_poly = gsw_t_freezing_poly(SA_chck_cast,p_chck_cast,0.5);

gsw_cv.pot_enthalpy_ice_freezing = gsw_pot_enthalpy_ice_freezing(SA_chck_cast,p_chck_cast);

gsw_cv.pot_enthalpy_ice_freezing_poly = gsw_pot_enthalpy_ice_freezing_poly(SA_chck_cast,p_chck_cast);

gsw_cv.SA_freezing_from_CT = gsw_SA_freezing_from_CT(gsw_cv.CT_freezing,p_chck_cast,0.5);

gsw_cv.SA_freezing_from_CT_poly = gsw_SA_freezing_from_CT_poly(gsw_cv.CT_freezing_poly,p_chck_cast,0.5);

gsw_cv.SA_freezing_from_t = gsw_SA_freezing_from_t(gsw_cv.t_freezing,p_chck_cast,0.5);

gsw_cv.SA_freezing_from_t_poly = gsw_SA_freezing_from_t_poly(gsw_cv.t_freezing_poly,p_chck_cast,0.5);

gsw_cv.pressure_freezing_CT = gsw_pressure_freezing_CT(gsw_cv.SA_Arctic,(gsw_cv.CT_Arctic - 1),0.5);

[gsw_cv.CTfreezing_SA, gsw_cv.CTfreezing_P] = gsw_CT_freezing_first_derivatives(SA_chck_cast,p_chck_cast,0.5);

[gsw_cv.CTfreezing_SA_poly, gsw_cv.CTfreezing_P_poly] = gsw_CT_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast,0.5);

[gsw_cv.tfreezing_SA, gsw_cv.tfreezing_P] = gsw_t_freezing_first_derivatives(SA_chck_cast,p_chck_cast,0.5);

[gsw_cv.tfreezing_SA_poly, gsw_cv.tfreezing_P_poly] = gsw_t_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast,0.5);

[gsw_cv.pot_enthalpy_ice_freezing_SA, gsw_cv.pot_enthalpy_ice_freezing_P] = gsw_pot_enthalpy_ice_freezing_first_derivatives(SA_chck_cast,p_chck_cast);

[gsw_cv.pot_enthalpy_ice_freezing_SA_poly, gsw_cv.pot_enthalpy_ice_freezing_P_poly] = gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast);

gsw_cv.latentheat_melting = gsw_latentheat_melting(SA_chck_cast,p_chck_cast);

%%  thermodynamic interaction between ice and seawater

gsw_cv.melting_ice_SA_CT_ratio = gsw_melting_ice_SA_CT_ratio(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.t_seaice);

gsw_cv.melting_ice_SA_CT_ratio_poly = gsw_melting_ice_SA_CT_ratio_poly(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.t_seaice);
 
gsw_cv.melting_ice_equilibrium_SA_CT_ratio = gsw_melting_ice_equilibrium_SA_CT_ratio(gsw_cv.SA_Arctic,gsw_cv.p_Arctic);

gsw_cv.melting_ice_equilibrium_SA_CT_ratio_poly = gsw_melting_ice_equilibrium_SA_CT_ratio_poly(gsw_cv.SA_Arctic,gsw_cv.p_Arctic);
                                           
[gsw_cv.melting_ice_into_seawater_SA_final, gsw_cv.melting_ice_into_seawater_CT_final] = gsw_melting_ice_into_seawater(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic+0.1,gsw_cv.p_Arctic,gsw_cv.w_seaice,gsw_cv.t_seaice);
 
[gsw_cv.ice_fraction_to_freeze_seawater_SA_freeze, gsw_cv.ice_fraction_to_freeze_seawater_CT_freeze, gsw_cv.ice_fraction_to_freeze_seawater_w_Ih]  = gsw_ice_fraction_to_freeze_seawater(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.t_seaice);
 
[gsw_cv.dSA_dCT_frazil, gsw_cv.dSA_dP_frazil, gsw_cv.dCT_dP_frazil] = gsw_frazil_ratios_adiabatic(gsw_cv.SA_Arctic,gsw_cv.p_Arctic,gsw_cv.w_seaice);

[gsw_cv.dSA_dCT_frazil_poly, gsw_cv.dSA_dP_frazil_poly, gsw_cv.dCT_dP_frazil_poly] = gsw_frazil_ratios_adiabatic_poly(gsw_cv.SA_Arctic,gsw_cv.p_Arctic,gsw_cv.w_seaice);

[gsw_cv.frazil_properties_potential_SA_final, gsw_cv.frazil_properties_potential_CT_final, gsw_cv.frazil_properties_potential_w_Ih_final] = gsw_frazil_properties_potential(gsw_cv.SA_bulk,gsw_cv.h_pot_bulk,gsw_cv.p_Arctic);

[gsw_cv.frazil_properties_potential_poly_SA_final, gsw_cv.frazil_properties_potential_poly_CT_final, gsw_cv.frazil_properties_potential_poly_w_Ih_final] = gsw_frazil_properties_potential_poly(gsw_cv.SA_bulk,gsw_cv.h_pot_bulk,gsw_cv.p_Arctic);

[gsw_cv.frazil_properties_SA_final, gsw_cv.frazil_properties_CT_final, gsw_cv.frazil_properties_w_Ih_final] = gsw_frazil_properties(gsw_cv.SA_bulk,gsw_cv.h_bulk,gsw_cv.p_Arctic);

%%  thermodynamic interaction between sea ice and seawater

gsw_cv.melting_seaice_SA_CT_ratio = gsw_melting_seaice_SA_CT_ratio(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.SA_seaice,gsw_cv.t_seaice);

gsw_cv.melting_seaice_SA_CT_ratio_poly = gsw_melting_seaice_SA_CT_ratio_poly(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.SA_seaice,gsw_cv.t_seaice);

gsw_cv.melting_seaice_equilibrium_SA_CT_ratio = gsw_melting_seaice_equilibrium_SA_CT_ratio(gsw_cv.SA_Arctic,gsw_cv.p_Arctic);   

gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_poly = gsw_melting_seaice_equilibrium_SA_CT_ratio_poly(gsw_cv.SA_Arctic,gsw_cv.p_Arctic);   
                                             
[gsw_cv.melting_seaice_into_seawater_SA_final, gsw_cv.melting_seaice_into_seawater_CT_final] = gsw_melting_seaice_into_seawater(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.w_seaice,gsw_cv.SA_seaice,gsw_cv.t_seaice);  

[gsw_cv.seaice_fraction_to_freeze_seawater_SA_freeze, gsw_cv.seaice_fraction_to_freeze_seawater_CT_freeze, gsw_cv.seaice_fraction_to_freeze_seawater_w_Ih]  = gsw_seaice_fraction_to_freeze_seawater(gsw_cv.SA_Arctic,gsw_cv.CT_Arctic,gsw_cv.p_Arctic,gsw_cv.SA_seaice,gsw_cv.t_seaice);  

%% themodynamic properties of ice Ih

gsw_cv.rho_ice = gsw_rho_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.alpha_wrt_t_ice = gsw_alpha_wrt_t_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.specvol_ice = gsw_specvol_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.pressure_coefficient_ice = gsw_pressure_coefficient_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.sound_speed_ice = gsw_sound_speed_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.kappa_ice = gsw_kappa_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.kappa_const_t_ice = gsw_kappa_const_t_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.internal_energy_ice = gsw_internal_energy_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.enthalpy_ice = gsw_enthalpy_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.entropy_ice = gsw_entropy_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);

gsw_cv.cp_ice = gsw_cp_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.chem_potential_water_ice = gsw_chem_potential_water_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.Helmholtz_energy_ice = gsw_Helmholtz_energy_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.adiabatic_lapse_rate_ice = gsw_adiabatic_lapse_rate_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.pt0_from_t_ice = gsw_pt0_from_t_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic);
 
gsw_cv.pt_from_t_ice = gsw_pt_from_t_ice(gsw_cv.t_seaice,gsw_cv.p_Arctic,pr);
 
gsw_cv.t_from_pt0_ice = gsw_t_from_pt0_ice(gsw_cv.pt0_from_t_ice,gsw_cv.p_Arctic);

gsw_cv.t_from_rho_ice = gsw_t_from_rho_ice(gsw_cv.rho_ice,gsw_cv.p_Arctic);

gsw_cv.pot_enthalpy_from_pt_ice = gsw_pot_enthalpy_from_pt_ice(gsw_cv.pt0_from_t_ice);
 
gsw_cv.pt_from_pot_enthalpy_ice = gsw_pt_from_pot_enthalpy_ice (gsw_cv.pot_enthalpy_from_pt_ice);
 
gsw_cv.pot_enthalpy_from_pt_ice_poly = gsw_pot_enthalpy_from_pt_ice_poly(gsw_cv.pt0_from_t_ice);
 
gsw_cv.pt_from_pot_enthalpy_ice_poly = gsw_pt_from_pot_enthalpy_ice_poly(gsw_cv.pot_enthalpy_from_pt_ice_poly);

gsw_cv.pot_enthalpy_from_specvol_ice = gsw_pot_enthalpy_from_specvol_ice(gsw_cv.specvol_ice,gsw_cv.p_Arctic);
 
gsw_cv.specvol_from_pot_enthalpy_ice = gsw_specvol_from_pot_enthalpy_ice (gsw_cv.pot_enthalpy_from_specvol_ice,gsw_cv.p_Arctic);
 
gsw_cv.pot_enthalpy_from_specvol_ice_poly = gsw_pot_enthalpy_from_specvol_ice_poly(gsw_cv.specvol_ice,gsw_cv.p_Arctic);
 
gsw_cv.specvol_from_pot_enthalpy_ice_poly = gsw_specvol_from_pot_enthalpy_ice_poly(gsw_cv.pot_enthalpy_from_specvol_ice_poly,gsw_cv.p_Arctic);


%% Isobaric isobaric evaporation enthalpy

gsw_cv.latentheat_evap_CT = gsw_latentheat_evap_CT(SA_chck_cast,CT_chck_cast);

gsw_cv.latentheat_evap_t = gsw_latentheat_evap_t(SA_chck_cast,t_chck_cast);

%% spiciness

gsw_cv.spiciness0 = gsw_spiciness0(SA_chck_cast,CT_chck_cast);

gsw_cv.spiciness1 = gsw_spiciness1(SA_chck_cast,CT_chck_cast);

gsw_cv.spiciness2 = gsw_spiciness2(SA_chck_cast,CT_chck_cast);

[gsw_cv.SA_spicsig0, gsw_cv.CT_spicsig0] = gsw_SA_CT_from_sigma0_spiciness0(gsw_cv.sigma0,gsw_cv.spiciness0);

[gsw_cv.SA_spicsig1, gsw_cv.CT_spicsig1] = gsw_SA_CT_from_sigma1_spiciness1(gsw_cv.sigma1,gsw_cv.spiciness1);

[gsw_cv.SA_spicsig2, gsw_cv.CT_spicsig2] = gsw_SA_CT_from_sigma2_spiciness2(gsw_cv.sigma2,gsw_cv.spiciness2);

%% planet earth properties 

gsw_cv.f = gsw_f(lat_chck_cast);
 
gsw_cv.grav = gsw_grav(lat_chck_cast,p_chck_cast);
 
gsw_cv.distance = gsw_distance(long_chck_cast,lat_chck_cast,p_chck_cast);  

%% TEOS-10 constants

gsw_cv.T0 = gsw_T0;

gsw_cv.P0 = gsw_P0;

gsw_cv.SSO = gsw_SSO;

gsw_cv.uPS = gsw_uPS;

gsw_cv.cp0 = gsw_cp0;

gsw_cv.C3515 = gsw_C3515;

gsw_cv.SonCl = gsw_SonCl;

gsw_cv.valence_factor = gsw_valence_factor;

gsw_cv.atomic_weight = gsw_atomic_weight;

%% dissolved gasses

gsw_cv.Arsol = gsw_Arsol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
 
gsw_cv.Arsol_SP_pt = gsw_Arsol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t);
 
gsw_cv.Hesol = gsw_Hesol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
 
gsw_cv.Hesol_SP_pt = gsw_Hesol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t); 
 
gsw_cv.Krsol = gsw_Krsol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);  
 
gsw_cv.Krsol_SP_pt = gsw_Krsol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t);   
 
gsw_cv.N2Osol = gsw_N2Osol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);   
 
gsw_cv.N2Osol_SP_pt = gsw_N2Osol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t); 
 
gsw_cv.N2sol = gsw_N2sol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast); 
 
gsw_cv.N2sol_SP_pt = gsw_N2sol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t);
 
gsw_cv.Nesol = gsw_Nesol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);   
 
gsw_cv.Nesol_SP_pt = gsw_Nesol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t);  
 
gsw_cv.O2sol = gsw_O2sol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);  
 
gsw_cv.O2sol_SP_pt = gsw_O2sol_SP_pt(SP_chck_cast,gsw_cv.pt0_from_t);  

%% density and enthalpy in terms of CT, derived from the exact Gibbs function

gsw_cv.specvol_CT_exact = gsw_specvol_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.alpha_CT_exact = gsw_alpha_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.beta_CT_exact = gsw_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.alpha_on_beta_CT_exact = gsw_alpha_on_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_vab_CT_exact, gsw_cv.alpha_vab_CT_exact, gsw_cv.beta_vab_CT_exact] = gsw_specvol_alpha_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
[gsw_cv.v_SA_CT_exact, gsw_cv.v_CT_CT_exact, gsw_cv.v_P_CT_exact] = gsw_specvol_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_SA_CT_exact, gsw_cv.v_SA_CT_CT_exact, gsw_cv.v_CT_CT_CT_exact, gsw_cv.v_SA_P_CT_exact, gsw_cv.v_CT_P_CT_exact] = gsw_specvol_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_wrt_h_CT_exact, gsw_cv.v_h_CT_exact] = gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.v_SA_SA_wrt_h_CT_exact, gsw_cv.v_SA_h_CT_exact, gsw_cv.v_h_h_CT_exact] = gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.specvol_anom_CT_exact = gsw_specvol_anom_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast,SA_ref,CT_ref);
 
gsw_cv.specvol_anom_standard_CT_exact = gsw_specvol_anom_standard_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
rho_CT_exact = gsw_rho_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
 gsw_cv.rho_CT_exact = rho_CT_exact;

[gsw_cv.rho_rab_CT_exact, gsw_cv.alpha_rab_CT_exact, gsw_cv.beta_rab_CT_exact] = gsw_rho_alpha_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_CT_exact, gsw_cv.rho_CT_CT_exact, gsw_cv.rho_P_CT_exact] = gsw_rho_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_SA_CT_exact, gsw_cv.rho_SA_CT_CT_exact, gsw_cv.rho_CT_CT_CT_exact, gsw_cv.rho_SA_P_CT_exact, gsw_cv.rho_CT_P_CT_exact] = gsw_rho_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_wrt_h_CT_exact, gsw_cv.rho_h_CT_exact] = gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.rho_SA_SA_wrt_h_CT_exact, gsw_cv.rho_SA_h_CT_exact, gsw_cv.rho_h_h_CT_exact] = gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.sigma0_CT_exact = gsw_sigma0_CT_exact(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma1_CT_exact = gsw_sigma1_CT_exact(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma2_CT_exact = gsw_sigma2_CT_exact(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma3_CT_exact = gsw_sigma3_CT_exact(SA_chck_cast,CT_chck_cast);
 
gsw_cv.sigma4_CT_exact = gsw_sigma4_CT_exact(SA_chck_cast,CT_chck_cast);

gsw_cv.cabbeling_CT_exact = gsw_cabbeling_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
gsw_cv.thermobaric_CT_exact = gsw_thermobaric_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.enthalpy_CT_exact = gsw_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.enthalpy_diff_CT_exact =  gsw_enthalpy_diff_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast_shallow,p_chck_cast_deep);
 
gsw_cv.dynamic_enthalpy_CT_exact = gsw_dynamic_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
 
[gsw_cv.h_SA_CT_exact, gsw_cv.h_CT_CT_exact] = gsw_enthalpy_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);  
 
[gsw_cv.h_SA_SA_CT_exact, gsw_cv.h_SA_CT_CT_exact, gsw_cv.h_CT_CT_CT_exact] = gsw_enthalpy_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast); 

gsw_cv.sound_speed_CT_exact = gsw_sound_speed_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.kappa_CT_exact = gsw_kappa_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.internal_energy_CT_exact = gsw_internal_energy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.u_SA_CT_exact, gsw_cv.u_CT_CT_exact, gsw_cv.u_P_CT_exact] = gsw_internal_energy_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

[gsw_cv.u_SA_SA_CT_exact, gsw_cv.u_SA_CT_CT_exact, gsw_cv.u_CT_CT_CT_exact, gsw_cv.u_SA_P_CT_exact, gsw_cv.u_CT_P_CT_exact] = gsw_internal_energy_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);

gsw_cv.CT_from_enthalpy_exact = gsw_CT_from_enthalpy_exact(SA_chck_cast,gsw_cv.enthalpy_CT_exact,p_chck_cast);

gsw_cv.SA_from_rho_CT_exact = gsw_SA_from_rho_CT_exact(rho_CT_exact,CT_chck_cast,p_chck_cast);

gsw_cv.CT_from_rho_exact = gsw_CT_from_rho_exact(rho_CT_exact,SA_chck_cast,p_chck_cast);

gsw_cv.CT_maxdensity_exact = gsw_CT_maxdensity_exact(SA_chck_cast,p_chck_cast);

%% Labrortory functions

rho_t_exact = gsw_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
gsw_cv.rho_t_exact = rho_t_exact;

gsw_cv.SA_from_rho_t_exact = gsw_SA_from_rho_t_exact(rho_t_exact,t_chck_cast,p_chck_cast);

gsw_cv.deltaSA_from_rho_t_exact = gsw_deltaSA_from_rho_t_exact(rho_t_exact,SP_chck_cast,t_chck_cast,p_chck_cast);

%% basic thermodynamic properties interms of in-situ t, derived from the exact Gibbs function

% rho_t_exact = gsw_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
% gsw_cv.rho_t_exact = rho_t_exact;

gsw_cv.pot_rho_t_exact = gsw_pot_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast,pr);

gsw_cv.sigma0_pt0_exact = gsw_sigma0_pt0_exact(SA_chck_cast,pt0);

gsw_cv.alpha_wrt_CT_t_exact = gsw_alpha_wrt_CT_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.alpha_wrt_pt_t_exact = gsw_alpha_wrt_pt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.alpha_wrt_t_exact = gsw_alpha_wrt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.beta_const_CT_t_exact = gsw_beta_const_CT_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.beta_const_pt_t_exact = gsw_beta_const_pt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.beta_const_t_exact = gsw_beta_const_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.specvol_t_exact = gsw_specvol_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast); 
 
gsw_cv.specvol_anom_standard_t_exact = gsw_specvol_anom_standard_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.sound_speed_t_exact = gsw_sound_speed_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.kappa_t_exact = gsw_kappa_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.kappa_const_t_exact = gsw_kappa_const_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.SA_from_rho_t_exact = gsw_SA_from_rho_t_exact(rho_t_exact,t_chck_cast,p_chck_cast);
 
gsw_cv.t_from_rho_exact = gsw_t_from_rho_exact(rho_t_exact,SA_chck_cast,p_chck_cast);

gsw_cv.t_maxdensity_exact = gsw_t_maxdensity_exact(SA_chck_cast,p_chck_cast);

gsw_cv.internal_energy_t_exact = gsw_internal_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.enthalpy_t_exact = gsw_enthalpy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.dynamic_enthalpy_t_exact = gsw_dynamic_enthalpy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

[gsw_cv.CT_SA_wrt_t, gsw_cv.CT_T_wrt_t, gsw_cv.CT_P_wrt_t] = gsw_CT_first_derivatives_wrt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);  

[gsw_cv.h_SA_wrt_t, gsw_cv.h_T_wrt_t, gsw_cv.h_P_wrt_t] = gsw_enthalpy_first_derivatives_wrt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);  

gsw_cv.cp_t_exact = gsw_cp_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.isochoric_heat_cap_t_exact = gsw_isochoric_heat_cap_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.chem_potential_relative_t_exact =  gsw_chem_potential_relative_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.chem_potential_water_t_exact =  gsw_chem_potential_water_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.chem_potential_salt_t_exact =  gsw_chem_potential_salt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.t_deriv_chem_potential_water_t_exact =  gsw_t_deriv_chem_potential_water_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.dilution_coefficient_t_exact = gsw_dilution_coefficient_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.Gibbs_energy_t_exact = gsw_Gibbs_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);

gsw_cv.Helmholtz_energy_t_exact = gsw_Helmholtz_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.osmotic_coefficient_t_exact = gsw_osmotic_coefficient_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 
gsw_cv.osmotic_pressure_t_exact = gsw_osmotic_pressure_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
 

%% Library

gsw_cv.Fdelta = gsw_Fdelta(p_chck_cast,long_chck_cast,lat_chck_cast);

for I = 1:45
    long_chck_cast_temp(I,:) = long_chck_cast(1,:);
    lat_chck_cast_temp(I,:) = lat_chck_cast(1,:);
end
[I] = find(~isnan(p_chck_cast));
gsw_cv.deltaSA_atlas = nan(45,3);
gsw_cv.deltaSA_atlas(I) = gsw_deltaSA_atlas(p_chck_cast(I),long_chck_cast_temp(I),lat_chck_cast_temp(I));



%save library\gsw_chck_vals_file.mat gsw_cv
save gsw_chck_vals_file_v3_06_12.mat gsw_cv


