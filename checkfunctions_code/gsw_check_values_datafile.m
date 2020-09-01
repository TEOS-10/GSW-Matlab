load gsw_data_v3_0

clear gsw_cv 

load ('gsw_chck_vals_errors_v3_06_12.mat','gsw_cv')

version_number = ('3.06.12');

save new_gsw_data_v3_06_12.mat -v6

% create a netcdf file.
%Open the file
mode = netcdf.getConstant('CLOBBER');
ncid = netcdf.create('.\gsw_data_v3_0.nc',mode);

varid = netcdf.getConstant('GLOBAL');

netcdf.putAtt(ncid,varid,'version_date',version_date);
netcdf.putAtt(ncid,varid,'version_number',version_number);
netcdf.putAtt(ncid,varid,'history','Created by Paul M. Barker on 6th July 2020 for the Gibbs-Seawater Oceanographic Toolbox');
netcdf.putAtt(ncid,varid,'title','Contains the global data set of Absolute Salinity Anomaly Ratio, the global data set of Absolute Salinity Anomaly Atlas, a reference cast and check values for all of the functions');

%Define the dimensions
dims_p_ref = netcdf.defDim(ncid,'nz',length(p_ref));
dims_lats_ref = netcdf.defDim(ncid,'ny',length(lats_ref));
dims_longs_ref = netcdf.defDim(ncid,'nx',length(longs_ref));

dims_p_ref_cast = netcdf.defDim(ncid,'nz_cast',length(p_ref_cast));
dims_lat_ref_cast = netcdf.defDim(ncid,'ny_cast',length(lat_ref_cast));
dims_long_ref_cast = netcdf.defDim(ncid,'nx_cast',length(long_ref_cast));

[m, n] = size(gsw_cv.SP_chck_cast);
dims_test_cast_m = netcdf.defDim(ncid,'test_cast_length',m);
dims_test_cast_n = netcdf.defDim(ncid,'test_cast_number',n);

[ma, na] = size(gsw_cv.SA_Arctic);
dims_Arctic_test_cast_m = netcdf.defDim(ncid,'Arctic_test_cast_length',ma);
dims_Arctic_test_cast_n = netcdf.defDim(ncid,'Arctic_test_cast_number',na);

dims_value_test_cast = netcdf.defDim(ncid,'value_test_cast',length(gsw_cv.pr));

[mm, nm] = size(gsw_cv.n2);
dims_test_cast_mm = netcdf.defDim(ncid,'test_cast_midpressure_length',mm);
dims_test_cast_nm = netcdf.defDim(ncid,'test_cast_midpressure_number',nm);

[ml, nl] = size(gsw_cv.distance);
dims_test_cast_ml = netcdf.defDim(ncid,'test_cast_midlocation_length',ml);
dims_test_cast_nl = netcdf.defDim(ncid,'test_cast_midlocation_number',nl);

[mi, ni] = size(gsw_cv.p_i);
dims_test_cast_mi = netcdf.defDim(ncid,'interp_test_cast_length',mi);
dims_test_cast_ni = netcdf.defDim(ncid,'interp_test_cast_number',ni);


%Define IDs for the dimension variables (pressure,latitude,longitude,...)
p_ref_ID = netcdf.defVar(ncid,'p_ref','double',dims_p_ref);
netcdf.putAtt(ncid,p_ref_ID,'standard_name','pressure_reference_atlas');
netcdf.putAtt(ncid,p_ref_ID,'units','dbar');

lats_ref_ID = netcdf.defVar(ncid,'lats_ref','double',dims_lats_ref);
netcdf.putAtt(ncid,lats_ref_ID,'standard_name','latitude_reference_atlas');
netcdf.putAtt(ncid,lats_ref_ID,'units','degrees_North');

longs_ref_ID = netcdf.defVar(ncid,'longs_ref','double',dims_longs_ref);
netcdf.putAtt(ncid,longs_ref_ID,'standard_name','longitude_reference_atlas');
netcdf.putAtt(ncid,longs_ref_ID,'units','degrees_East');

p_ref_cast_ID = netcdf.defVar(ncid,'p_ref_cast','double',dims_p_ref_cast);
netcdf.putAtt(ncid,p_ref_cast_ID,'standard_name','pressure_reference_cast');
netcdf.putAtt(ncid,p_ref_cast_ID,'units','dbar');

lat_ref_cast_ID = netcdf.defVar(ncid,'lat_ref_cast','double',dims_lat_ref_cast);
netcdf.putAtt(ncid,lat_ref_cast_ID,'standard_name','latitude_reference_cast');
netcdf.putAtt(ncid,lat_ref_cast_ID,'units','degrees_North');

long_ref_cast_ID = netcdf.defVar(ncid,'long_ref_cast','double',dims_long_ref_cast);
netcdf.putAtt(ncid,long_ref_cast_ID,'standard_name','longitude_reference_cast');
netcdf.putAtt(ncid,long_ref_cast_ID,'units','degrees_East');

%Define the main variables (SAAR,SR,delta_SA,ocean,ndepth)
SAAR_ref_ID = netcdf.defVar(ncid,'SAAR_ref','double',[dims_p_ref dims_lats_ref dims_longs_ref]);
netcdf.putAtt(ncid,SAAR_ref_ID,'standard_name','Absolute_Salinity_Anomaly_Ratio');
netcdf.putAtt(ncid,SAAR_ref_ID,'units','unitless');
netcdf.putAtt(ncid,SAAR_ref_ID,'non_ocean_value',-9.e-99);

SR_ref_ID = netcdf.defVar(ncid,'SA_ref','double',[dims_p_ref dims_lats_ref dims_longs_ref]);
netcdf.putAtt(ncid,SR_ref_ID,'standard_name','Reference_Salinity_atlas_value');
netcdf.putAtt(ncid,SR_ref_ID,'units','g/kg');
netcdf.putAtt(ncid,SR_ref_ID,'non_ocean_value',-9.e-99);

deltaSA_ref_ID = netcdf.defVar(ncid,'deltaSA_ref','double',[dims_p_ref dims_lats_ref dims_longs_ref]);
netcdf.putAtt(ncid,deltaSA_ref_ID,'standard_name','Absolute_Salinity_Anomaly_atlas_value');
netcdf.putAtt(ncid,deltaSA_ref_ID,'units','g/kg');
netcdf.putAtt(ncid,deltaSA_ref_ID,'non_ocean_value',-9.e-99);

ocean_ref_ID = netcdf.defVar(ncid,'ocean_ref','double',[dims_lats_ref dims_longs_ref]);
netcdf.putAtt(ncid,ocean_ref_ID,'standard_name','ocean_atlas_value_in_atlas_at_that_grid');
netcdf.putAtt(ncid,ocean_ref_ID,'units','unitless');
netcdf.putAtt(ncid,ocean_ref_ID,'non_ocean_value',-9.e-99);

ndepth_ref_ID = netcdf.defVar(ncid,'ndepth_ref','double',[dims_lats_ref dims_longs_ref]);
netcdf.putAtt(ncid,ndepth_ref_ID,'standard_name','number of levels_in_atlas_at_that_grid');

SA_ref_cast_ID = netcdf.defVar(ncid,'SA_ref_cast','double',dims_p_ref_cast);
netcdf.putAtt(ncid,SA_ref_cast_ID,'standard_name','Absolute_Salinity_reference_cast');
netcdf.putAtt(ncid,SA_ref_cast_ID,'units','g/kg');

CT_ref_cast_ID = netcdf.defVar(ncid,'CT_ref_cast','double',dims_p_ref_cast);
netcdf.putAtt(ncid,CT_ref_cast_ID,'standard_name','Conservative_Temperature_reference_cast');
netcdf.putAtt(ncid,CT_ref_cast_ID,'units','degrees_Celcius');

gamma_n_ref_cast_ID = netcdf.defVar(ncid,'gamma_n_ref_cast','double',dims_p_ref_cast);
netcdf.putAtt(ncid,gamma_n_ref_cast_ID,'standard_name','gamma_n_reference_cast');
netcdf.putAtt(ncid,gamma_n_ref_cast_ID,'units','g/kg');

sigma_2_ref_cast_ID = netcdf.defVar(ncid,'sigma_2_ref_cast','double',dims_p_ref_cast);
netcdf.putAtt(ncid,sigma_2_ref_cast_ID,'standard_name','sigma_2_reference_cast');
netcdf.putAtt(ncid,sigma_2_ref_cast_ID,'units','g/m^3');

SA_Arctic_ID = netcdf.defVar(ncid,'SA_Arctic','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

t_Arctic_ID = netcdf.defVar(ncid,'t_Arctic','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

p_Arctic_ID = netcdf.defVar(ncid,'p_Arctic','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

latitude_Arctic_ID = netcdf.defVar(ncid,'latitude_Arctic','double',dims_Arctic_test_cast_n);

longitude_Arctic_ID = netcdf.defVar(ncid,'longitude_Arctic','double',dims_Arctic_test_cast_n);

SA_seaice_ID = netcdf.defVar(ncid,'SA_seaice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

t_seaice_ID = netcdf.defVar(ncid,'t_seaice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

w_seaice_ID = netcdf.defVar(ncid,'w_seaice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

CT_Arctic_ID = netcdf.defVar(ncid,'CT_Arctic','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

t_ice_ID = netcdf.defVar(ncid,'t_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

w_ice_ID = netcdf.defVar(ncid,'w_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

SA_bulk_ID = netcdf.defVar(ncid,'SA_bulk','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

h_bulk_ID = netcdf.defVar(ncid,'h_bulk','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

h_pot_bulk_ID = netcdf.defVar(ncid,'h_pot_bulk','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);

CT_chck_cast_ID = netcdf.defVar(ncid,'CT_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

Rt_chck_cast_ID = netcdf.defVar(ncid,'Rt_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

SA_chck_cast_ID = netcdf.defVar(ncid,'SA_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

SK_chck_cast_ID = netcdf.defVar(ncid,'SK_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

SP_chck_cast_ID = netcdf.defVar(ncid,'SP_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

t_chck_cast_ID = netcdf.defVar(ncid,'t_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

p_chck_cast_ID = netcdf.defVar(ncid,'p_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

lat_chck_cast_ID = netcdf.defVar(ncid,'lat_chck_cast','double',dims_test_cast_n);

long_chck_cast_ID = netcdf.defVar(ncid,'long_chck_cast','double',dims_test_cast_n);

pr_ID = netcdf.defVar(ncid,'pr','double',dims_value_test_cast);

pr_05_ID = netcdf.defVar(ncid,'pr_05','double',dims_value_test_cast);

p_chck_cast_shallow_ID = netcdf.defVar(ncid,'p_chck_cast_shallow','double',[dims_test_cast_m dims_test_cast_n]);

p_chck_cast_deep_ID = netcdf.defVar(ncid,'p_chck_cast_deep','double',[dims_test_cast_m dims_test_cast_n]);

delta_p_chck_cast_ID = netcdf.defVar(ncid,'delta_p_chck_cast','double',[dims_test_cast_m dims_test_cast_n]);

Neutral_Density_ID = netcdf.defVar(ncid,'Neutral_Density','double', dims_test_cast_n);

p_Neutral_Density_ID = netcdf.defVar(ncid,'p_Neutral_Density','double', dims_test_cast_n);

%=============================================================================

%% Practical Salinity (SP):- PSS-78

C_from_SP_ID = netcdf.defVar(ncid,'C_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,C_from_SP_ID,'computation_accuracy',gsw_cv.C_from_SP_ca);

SP_from_C_ID = netcdf.defVar(ncid,'SP_from_C','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_C_ID,'computation_accuracy',gsw_cv.SP_from_C_ca);

R_from_SP_ID = netcdf.defVar(ncid,'R_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,R_from_SP_ID,'computation_accuracy',gsw_cv.R_from_SP_ca);

SP_from_R_ID = netcdf.defVar(ncid,'SP_from_R','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_R_ID,'computation_accuracy',gsw_cv.SP_from_R_ca);

SP_salinometer_ID = netcdf.defVar(ncid,'SP_salinometer','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_salinometer_ID,'computation_accuracy',gsw_cv.SP_salinometer_ca);

SP_from_SK_ID = netcdf.defVar(ncid,'SP_from_SK','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_SK_ID,'computation_accuracy',gsw_cv.SP_from_SK_ca);

%% Absolute Salinity (SA), Preformed Salinity (Sstar) and Conservative Temperature (CT)

SA_from_SP_ID = netcdf.defVar(ncid,'SA_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_SP_ID,'computation_accuracy',gsw_cv.SA_from_SP_ca);

Sstar_from_SP_ID = netcdf.defVar(ncid,'Sstar_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Sstar_from_SP_ID,'computation_accuracy',gsw_cv.Sstar_from_SP_ca);

CT_from_t_ID = netcdf.defVar(ncid,'CT_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_t_ID,'computation_accuracy',gsw_cv.CT_from_t_ca);

%% other conversions between temperatures, salinities, pressure and height

deltaSA_from_SP_ID = netcdf.defVar(ncid,'deltaSA_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,deltaSA_from_SP_ID,'computation_accuracy',gsw_cv.deltaSA_from_SP_ca);

SA_SA_Sstar_from_SP_ID = netcdf.defVar(ncid,'SA_SA_Sstar_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_SA_Sstar_from_SP_ID,'computation_accuracy',gsw_cv.SA_SA_Sstar_from_SP_ca);

Sstar_SA_Sstar_from_SP_ID = netcdf.defVar(ncid,'Sstar_SA_Sstar_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Sstar_SA_Sstar_from_SP_ID,'computation_accuracy',gsw_cv.Sstar_SA_Sstar_from_SP_ca);

SR_from_SP_ID = netcdf.defVar(ncid,'SR_from_SP','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SR_from_SP_ID,'computation_accuracy',gsw_cv.SR_from_SP_ca);

SP_from_SR_ID = netcdf.defVar(ncid,'SP_from_SR','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_SR_ID,'computation_accuracy',gsw_cv.SP_from_SR_ca);

SP_from_SA_ID = netcdf.defVar(ncid,'SP_from_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_SA_ID,'computation_accuracy',gsw_cv.SP_from_SA_ca);

Sstar_from_SA_ID = netcdf.defVar(ncid,'Sstar_from_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Sstar_from_SA_ID,'computation_accuracy',gsw_cv.Sstar_from_SA_ca);

SA_from_Sstar_ID = netcdf.defVar(ncid,'SA_from_Sstar','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_Sstar_ID,'computation_accuracy',gsw_cv.SA_from_Sstar_ca);

SP_from_Sstar_ID = netcdf.defVar(ncid,'SP_from_Sstar','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SP_from_Sstar_ID,'computation_accuracy',gsw_cv.SP_from_Sstar_ca);

t_from_CT_ID = netcdf.defVar(ncid,'t_from_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_from_CT_ID,'computation_accuracy',gsw_cv.t_from_CT_ca);

pt_from_CT_ID = netcdf.defVar(ncid,'pt_from_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_from_CT_ID,'computation_accuracy',gsw_cv.pt_from_CT_ca);

CT_from_pt_ID = netcdf.defVar(ncid,'CT_from_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_pt_ID,'computation_accuracy',gsw_cv.CT_from_pt_ca);

pot_enthalpy_from_pt_ID = netcdf.defVar(ncid,'pot_enthalpy_from_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_from_pt_ID,'computation_accuracy',gsw_cv.pot_enthalpy_from_pt_ca);

pt0_from_t_ID = netcdf.defVar(ncid,'pt0_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt0_from_t_ID,'computation_accuracy',gsw_cv.pt0_from_t_ca);

pt_from_t_ID = netcdf.defVar(ncid,'pt_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_from_t_ID,'computation_accuracy',gsw_cv.pt_from_t_ca);

t90_from_t68_ID = netcdf.defVar(ncid,'t90_from_t68','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t90_from_t68_ID,'computation_accuracy',gsw_cv.t90_from_t68_ca);

t90_from_t48_ID = netcdf.defVar(ncid,'t90_from_t48','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t90_from_t48_ID,'computation_accuracy',gsw_cv.t90_from_t48_ca);

z_from_p_ID = netcdf.defVar(ncid,'z_from_p','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,z_from_p_ID,'computation_accuracy',gsw_cv.z_from_p_ca);

p_from_z_ID = netcdf.defVar(ncid,'p_from_z','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,p_from_z_ID,'computation_accuracy',gsw_cv.p_from_z_ca);

depth_from_z_ID = netcdf.defVar(ncid,'depth_from_z','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,depth_from_z_ID,'computation_accuracy',gsw_cv.depth_from_z_ca);

z_from_depth_ID = netcdf.defVar(ncid,'z_from_depth','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,z_from_depth_ID,'computation_accuracy',gsw_cv.z_from_depth_ca);

Abs_Pressure_from_p_ID = netcdf.defVar(ncid,'Abs_Pressure_from_p','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Abs_Pressure_from_p_ID,'computation_accuracy',gsw_cv.Abs_Pressure_from_p_ca);

p_from_Abs_Pressure_ID = netcdf.defVar(ncid,'p_from_Abs_Pressure','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,p_from_Abs_Pressure_ID,'computation_accuracy',gsw_cv.p_from_Abs_Pressure_ca);

entropy_from_CT_ID = netcdf.defVar(ncid,'entropy_from_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,entropy_from_CT_ID,'computation_accuracy',gsw_cv.entropy_from_CT_ca);

CT_from_entropy_ID = netcdf.defVar(ncid,'CT_from_entropy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_entropy_ID,'computation_accuracy',gsw_cv.CT_from_entropy_ca);

entropy_from_pt_ID = netcdf.defVar(ncid,'entropy_from_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,entropy_from_pt_ID,'computation_accuracy',gsw_cv.entropy_from_pt_ca);

pt_from_entropy_ID = netcdf.defVar(ncid,'pt_from_entropy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_from_entropy_ID,'computation_accuracy',gsw_cv.pt_from_entropy_ca);

entropy_from_t_ID = netcdf.defVar(ncid,'entropy_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,entropy_from_t_ID,'computation_accuracy',gsw_cv.entropy_from_t_ca);

t_from_entropy_ID = netcdf.defVar(ncid,'t_from_entropy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_from_entropy_ID,'computation_accuracy',gsw_cv.t_from_entropy_ca);

adiabatic_lapse_rate_from_CT_ID = netcdf.defVar(ncid,'adiabatic_lapse_rate_from_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,adiabatic_lapse_rate_from_CT_ID,'computation_accuracy',gsw_cv.adiabatic_lapse_rate_from_CT_ca);

adiabatic_lapse_rate_from_t_ID = netcdf.defVar(ncid,'adiabatic_lapse_rate_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,adiabatic_lapse_rate_from_t_ID,'computation_accuracy',gsw_cv.adiabatic_lapse_rate_from_t_ca);

molality_from_SA_ID = netcdf.defVar(ncid,'molality_from_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,molality_from_SA_ID,'computation_accuracy',gsw_cv.molality_from_SA_ca);

ionic_strength_from_SA_ID = netcdf.defVar(ncid,'ionic_strength_from_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,ionic_strength_from_SA_ID,'computation_accuracy',gsw_cv.ionic_strength_from_SA_ca);

%% specific volume, density and enthalpy

specvol_ID = netcdf.defVar(ncid,'specvol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_ID,'computation_accuracy',gsw_cv.specvol_ca);

alpha_ID = netcdf.defVar(ncid,'alpha','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_ID,'computation_accuracy',gsw_cv.alpha_ca);

beta_ID = netcdf.defVar(ncid,'beta','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_ID,'computation_accuracy',gsw_cv.beta_ca);

alpha_on_beta_ID = netcdf.defVar(ncid,'alpha_on_beta','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_on_beta_ID,'computation_accuracy',gsw_cv.alpha_on_beta_ca);

v_vab_ID = netcdf.defVar(ncid,'v_vab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_vab_ID,'computation_accuracy',gsw_cv.v_vab_ca);

alpha_vab_ID = netcdf.defVar(ncid,'alpha_vab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_vab_ID,'computation_accuracy',gsw_cv.alpha_vab_ca);

beta_vab_ID = netcdf.defVar(ncid,'beta_vab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_vab_ID,'computation_accuracy',gsw_cv.beta_vab_ca);

v_SA_ID = netcdf.defVar(ncid,'v_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_ID,'computation_accuracy',gsw_cv.v_SA_ca);

v_CT_ID = netcdf.defVar(ncid,'v_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_ID,'computation_accuracy',gsw_cv.v_CT_ca);

v_P_ID = netcdf.defVar(ncid,'v_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_P_ID,'computation_accuracy',gsw_cv.v_P_ca);

v_SA_SA_ID = netcdf.defVar(ncid,'v_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_SA_ID,'computation_accuracy',gsw_cv.v_SA_SA_ca);

v_SA_CT_ID = netcdf.defVar(ncid,'v_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_CT_ID,'computation_accuracy',gsw_cv.v_SA_CT_ca);

v_CT_CT_ID = netcdf.defVar(ncid,'v_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_CT_ID,'computation_accuracy',gsw_cv.v_CT_CT_ca);

v_SA_P_ID = netcdf.defVar(ncid,'v_SA_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_P_ID,'computation_accuracy',gsw_cv.v_SA_P_ca);

v_CT_P_ID = netcdf.defVar(ncid,'v_CT_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_P_ID,'computation_accuracy',gsw_cv.v_CT_P_ca);

v_SA_wrt_h_ID = netcdf.defVar(ncid,'v_SA_wrt_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_wrt_h_ID,'computation_accuracy',gsw_cv.v_SA_wrt_h_ca);

v_h_ID = netcdf.defVar(ncid,'v_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_h_ID,'computation_accuracy',gsw_cv.v_h_ca);

v_SA_SA_wrt_h_ID = netcdf.defVar(ncid,'v_SA_SA_wrt_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_SA_wrt_h_ID,'computation_accuracy',gsw_cv.v_SA_SA_wrt_h_ca);

v_SA_h_ID = netcdf.defVar(ncid,'v_SA_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_h_ID,'computation_accuracy',gsw_cv.v_SA_h_ca);

v_h_h_ID = netcdf.defVar(ncid,'v_h_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_h_h_ID,'computation_accuracy',gsw_cv.v_h_h_ca);

specvol_anom_ID = netcdf.defVar(ncid,'specvol_anom','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_anom_ID,'computation_accuracy',gsw_cv.specvol_anom_ca);

specvol_anom_standard_ID = netcdf.defVar(ncid,'specvol_anom_standard','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_anom_standard_ID,'computation_accuracy',gsw_cv.specvol_anom_standard_ca);

rho_ID = netcdf.defVar(ncid,'rho','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_ID,'computation_accuracy',gsw_cv.rho_ca);

rho_rab_ID = netcdf.defVar(ncid,'rho_rab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_rab_ID,'computation_accuracy',gsw_cv.rho_rab_ca);

alpha_rab_ID = netcdf.defVar(ncid,'alpha_rab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_rab_ID,'computation_accuracy',gsw_cv.alpha_rab_ca);

beta_rab_ID = netcdf.defVar(ncid,'beta_rab','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_rab_ID,'computation_accuracy',gsw_cv.beta_rab_ca);

rho_SA_ID = netcdf.defVar(ncid,'rho_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_ID,'computation_accuracy',gsw_cv.rho_SA_ca);

rho_CT_ID = netcdf.defVar(ncid,'rho_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_ID,'computation_accuracy',gsw_cv.rho_CT_ca);

rho_P_ID = netcdf.defVar(ncid,'rho_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_P_ID,'computation_accuracy',gsw_cv.rho_P_ca);

rho_SA_SA_ID = netcdf.defVar(ncid,'rho_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_SA_ID,'computation_accuracy',gsw_cv.rho_SA_SA_ca);

rho_SA_CT_ID = netcdf.defVar(ncid,'rho_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_CT_ID,'computation_accuracy',gsw_cv.rho_SA_CT_ca);

rho_CT_CT_ID = netcdf.defVar(ncid,'rho_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_CT_ID,'computation_accuracy',gsw_cv.rho_CT_CT_ca);

rho_SA_P_ID = netcdf.defVar(ncid,'rho_SA_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_P_ID,'computation_accuracy',gsw_cv.rho_SA_P_ca);

rho_CT_P_ID = netcdf.defVar(ncid,'rho_CT_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_P_ID,'computation_accuracy',gsw_cv.rho_CT_P_ca);

rho_SA_wrt_h_ID = netcdf.defVar(ncid,'rho_SA_wrt_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_wrt_h_ID,'computation_accuracy',gsw_cv.rho_SA_wrt_h_ca);

rho_h_ID = netcdf.defVar(ncid,'rho_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_h_ID,'computation_accuracy',gsw_cv.rho_h_ca);

rho_SA_SA_wrt_h_ID = netcdf.defVar(ncid,'rho_SA_SA_wrt_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_SA_wrt_h_ID,'computation_accuracy',gsw_cv.rho_SA_SA_wrt_h_ca);

rho_SA_h_ID = netcdf.defVar(ncid,'rho_SA_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_h_ID,'computation_accuracy',gsw_cv.rho_SA_h_ca);

rho_h_h_ID = netcdf.defVar(ncid,'rho_h_h','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_h_h_ID,'computation_accuracy',gsw_cv.rho_h_h_ca);

sigma0_ID = netcdf.defVar(ncid,'sigma0','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma0_ID,'computation_accuracy',gsw_cv.sigma0_ca);

sigma1_ID = netcdf.defVar(ncid,'sigma1','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma1_ID,'computation_accuracy',gsw_cv.sigma1_ca);

sigma2_ID = netcdf.defVar(ncid,'sigma2','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma2_ID,'computation_accuracy',gsw_cv.sigma2_ca);

sigma3_ID = netcdf.defVar(ncid,'sigma3','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma3_ID,'computation_accuracy',gsw_cv.sigma3_ca);

sigma4_ID = netcdf.defVar(ncid,'sigma4','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma4_ID,'computation_accuracy',gsw_cv.sigma4_ca);

cabbeling_ID = netcdf.defVar(ncid,'cabbeling','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,cabbeling_ID,'computation_accuracy',gsw_cv.cabbeling_ca);

thermobaric_ID = netcdf.defVar(ncid,'thermobaric','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,thermobaric_ID,'computation_accuracy',gsw_cv.thermobaric_ca);

enthalpy_ID = netcdf.defVar(ncid,'enthalpy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_ID,'computation_accuracy',gsw_cv.enthalpy_ca);

enthalpy_diff_ID = netcdf.defVar(ncid,'enthalpy_diff','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_diff_ID,'computation_accuracy',gsw_cv.enthalpy_diff_ca);

dynamic_enthalpy_ID = netcdf.defVar(ncid,'dynamic_enthalpy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,dynamic_enthalpy_ID,'computation_accuracy',gsw_cv.dynamic_enthalpy_ca);

h_SA_ID = netcdf.defVar(ncid,'h_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_ID,'computation_accuracy',gsw_cv.h_SA_ca);

h_CT_ID = netcdf.defVar(ncid,'h_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_CT_ID,'computation_accuracy',gsw_cv.h_CT_ca);

h_SA_SA_ID = netcdf.defVar(ncid,'h_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_SA_ID,'computation_accuracy',gsw_cv.h_SA_SA_ca);

h_SA_CT_ID = netcdf.defVar(ncid,'h_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_CT_ID,'computation_accuracy',gsw_cv.h_SA_CT_ca);

h_CT_CT_ID = netcdf.defVar(ncid,'h_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_CT_CT_ID,'computation_accuracy',gsw_cv.h_CT_CT_ca);

sound_speed_ID = netcdf.defVar(ncid,'sound_speed','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sound_speed_ID,'computation_accuracy',gsw_cv.sound_speed_ca);

kappa_ID = netcdf.defVar(ncid,'kappa','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,kappa_ID,'computation_accuracy',gsw_cv.kappa_ca);

internal_energy_ID = netcdf.defVar(ncid,'internal_energy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,internal_energy_ID,'computation_accuracy',gsw_cv.internal_energy_ca);

u_SA_ID = netcdf.defVar(ncid,'u_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_ID,'computation_accuracy',gsw_cv.u_SA_ca);

u_CT_ID = netcdf.defVar(ncid,'u_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_ID,'computation_accuracy',gsw_cv.u_CT_ca);

u_P_ID = netcdf.defVar(ncid,'u_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_P_ID,'computation_accuracy',gsw_cv.u_P_ca);

u_SA_SA_ID = netcdf.defVar(ncid,'u_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_SA_ID,'computation_accuracy',gsw_cv.u_SA_SA_ca);

u_SA_CT_ID = netcdf.defVar(ncid,'u_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_CT_ID,'computation_accuracy',gsw_cv.u_SA_CT_ca);

u_CT_CT_ID = netcdf.defVar(ncid,'u_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_CT_ID,'computation_accuracy',gsw_cv.u_CT_CT_ca);

u_SA_P_ID = netcdf.defVar(ncid,'u_SA_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_P_ID,'computation_accuracy',gsw_cv.u_SA_P_ca);

u_CT_P_ID = netcdf.defVar(ncid,'u_CT_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_P_ID,'computation_accuracy',gsw_cv.u_CT_P_ca);

CT_from_enthalpy_ID = netcdf.defVar(ncid,'CT_from_enthalpy','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_enthalpy_ID,'computation_accuracy',gsw_cv.CT_from_enthalpy_ca);

SA_from_rho_ID = netcdf.defVar(ncid,'SA_from_rho','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_rho_ID,'computation_accuracy',gsw_cv.SA_from_rho_ca);

CT_from_rho_ID = netcdf.defVar(ncid,'CT_from_rho','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_rho_ID,'computation_accuracy',gsw_cv.CT_from_rho_ca);

CT_maxdensity_ID = netcdf.defVar(ncid,'CT_maxdensity','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_maxdensity_ID,'computation_accuracy',gsw_cv.CT_maxdensity_ca);

%% vertical stability and interpolation

Tu_ID = netcdf.defVar(ncid,'Tu','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,Tu_ID,'computation_accuracy',gsw_cv.Tu_ca);

Rsubrho_ID = netcdf.defVar(ncid,'Rsubrho','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,Rsubrho_ID,'computation_accuracy',gsw_cv.Rsubrho_ca);

p_mid_TuRsr_ID = netcdf.defVar(ncid,'p_mid_TuRsr','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,p_mid_TuRsr_ID,'computation_accuracy',gsw_cv.p_mid_TuRsr_ca);

n2_ID = netcdf.defVar(ncid,'n2','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2_ID,'computation_accuracy',gsw_cv.n2_ca);

p_mid_n2_ID = netcdf.defVar(ncid,'p_mid_n2','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,p_mid_n2_ID,'computation_accuracy',gsw_cv.p_mid_n2_ca);

n2min_ID = netcdf.defVar(ncid,'n2min','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_ID,'computation_accuracy',gsw_cv.n2min_ca);

n2min_pmid_ID = netcdf.defVar(ncid,'n2min_pmid','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_pmid_ID,'computation_accuracy',gsw_cv.n2min_pmid_ca);

n2min_specvol_ID = netcdf.defVar(ncid,'n2min_specvol','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_specvol_ID,'computation_accuracy',gsw_cv.n2min_specvol_ca);

n2min_alpha_ID = netcdf.defVar(ncid,'n2min_alpha','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_alpha_ID,'computation_accuracy',gsw_cv.n2min_alpha_ca);

n2min_beta_ID = netcdf.defVar(ncid,'n2min_beta','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_beta_ID,'computation_accuracy',gsw_cv.n2min_beta_ca);

n2min_dsa_ID = netcdf.defVar(ncid,'n2min_dsa','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_dsa_ID,'computation_accuracy',gsw_cv.n2min_dsa_ca);

n2min_dct_ID = netcdf.defVar(ncid,'n2min_dct','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_dct_ID,'computation_accuracy',gsw_cv.n2min_dct_ca);

n2min_dp_ID = netcdf.defVar(ncid,'n2min_dp','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,n2min_dp_ID,'computation_accuracy',gsw_cv.n2min_dp_ca);

mlp_ID = netcdf.defVar(ncid,'mlp','double',[dims_test_cast_n]);
netcdf.putAtt(ncid,mlp_ID,'computation_accuracy',gsw_cv.mlp_ca);

n2_lowerlimit_ID = netcdf.defVar(ncid,'n2_lowerlimit','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,n2_lowerlimit_ID,'computation_accuracy',gsw_cv.n2_lowerlimit_ca);


SAi_SACTinterp_ID = netcdf.defVar(ncid,'SAi_SACTinterp','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,SAi_SACTinterp_ID,'computation_accuracy',gsw_cv.SAi_SACTinterp_ca);

CTi_SACTinterp_ID = netcdf.defVar(ncid,'CTi_SACTinterp','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,CTi_SACTinterp_ID,'computation_accuracy',gsw_cv.CTi_SACTinterp_ca);

ti_tinterp_ID = netcdf.defVar(ncid,'ti_tinterp','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,ti_tinterp_ID,'computation_accuracy',gsw_cv.ti_tinterp_ca);

traceri_tracerCTinterp_ID = netcdf.defVar(ncid,'traceri_tracerCTinterp','double',[dims_test_cast_mi dims_test_cast_ni]);
netcdf.putAtt(ncid,traceri_tracerCTinterp_ID,'computation_accuracy',gsw_cv.traceri_tracerCTinterp_ca);

CTi_tracerCTinterp_ID = netcdf.defVar(ncid,'CTi_tracerCTinterp','double',[dims_test_cast_mi dims_test_cast_ni]);
netcdf.putAtt(ncid,CTi_tracerCTinterp_ID,'computation_accuracy',gsw_cv.CTi_tracerCTinterp_ca);

traceri_tracerinterp_ID = netcdf.defVar(ncid,'traceri_tracerinterp','double',[dims_test_cast_mi dims_test_cast_ni]);
netcdf.putAtt(ncid,traceri_tracerinterp_ID,'computation_accuracy',gsw_cv.traceri_tracerinterp_ca);


IPVfN2_ID = netcdf.defVar(ncid,'IPVfN2','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,IPVfN2_ID,'computation_accuracy',gsw_cv.IPVfN2_ca);

p_mid_IPVfN2_ID = netcdf.defVar(ncid,'p_mid_IPVfN2','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,p_mid_IPVfN2_ID,'computation_accuracy',gsw_cv.p_mid_IPVfN2_ca);

%% geostrophic streamfunctions and acoustic travel time

geo_strf_dyn_height_ID = netcdf.defVar(ncid,'geo_strf_dyn_height','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_dyn_height_ID,'computation_accuracy',gsw_cv.geo_strf_dyn_height_ca);

geo_strf_dyn_height_pc_ID = netcdf.defVar(ncid,'geo_strf_dyn_height_pc','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_dyn_height_pc_ID,'computation_accuracy',gsw_cv.geo_strf_dyn_height_pc_ca);

geo_strf_dyn_height_pc_p_mid_ID = netcdf.defVar(ncid,'geo_strf_dyn_height_pc_p_mid','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_dyn_height_pc_p_mid_ID,'computation_accuracy',gsw_cv.geo_strf_dyn_height_pc_p_mid_ca);

geo_strf_isopycnal_ID = netcdf.defVar(ncid,'geo_strf_isopycnal','double',dims_test_cast_n);
netcdf.putAtt(ncid,geo_strf_isopycnal_ID,'computation_accuracy',gsw_cv.geo_strf_isopycnal_ca);

geo_strf_isopycnal_pc_ID = netcdf.defVar(ncid,'geo_strf_isopycnal_pc','double',dims_test_cast_n);
netcdf.putAtt(ncid,geo_strf_isopycnal_pc_ID,'computation_accuracy',gsw_cv.geo_strf_isopycnal_pc_ca);

geo_strf_isopycnal_pc_p_mid_ID = netcdf.defVar(ncid,'geo_strf_isopycnal_pc_p_mid','double', dims_test_cast_n);
netcdf.putAtt(ncid,geo_strf_isopycnal_pc_p_mid_ID,'computation_accuracy',gsw_cv.geo_strf_isopycnal_pc_p_mid_ca);

geo_strf_Montgomery_ID = netcdf.defVar(ncid,'geo_strf_Montgomery','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_Montgomery_ID,'computation_accuracy',gsw_cv.geo_strf_Montgomery_ca);

geo_strf_Cunningham_ID = netcdf.defVar(ncid,'geo_strf_Cunningham','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_Cunningham_ID,'computation_accuracy',gsw_cv.geo_strf_Cunningham_ca);

geo_strf_steric_height_ID = netcdf.defVar(ncid,'geo_strf_steric_height','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_steric_height_ID,'computation_accuracy',gsw_cv.geo_strf_steric_height_ca);

geo_strf_PISH_ID = netcdf.defVar(ncid,'geo_strf_PISH','double',[dims_test_cast_n]);
netcdf.putAtt(ncid,geo_strf_PISH_ID,'computation_accuracy',gsw_cv.geo_strf_PISH_ca);

travel_time_ID = netcdf.defVar(ncid,'travel_time','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,travel_time_ID,'computation_accuracy',gsw_cv.travel_time_ca);

%% Geostrophic velocity

geo_strf_velocity_ID = netcdf.defVar(ncid,'geo_strf_velocity','double',[dims_test_cast_ml dims_test_cast_nl]);
netcdf.putAtt(ncid,geo_strf_velocity_ID,'computation_accuracy',gsw_cv.geo_strf_velocity_ca);

geo_strf_velocity_mid_lat_ID = netcdf.defVar(ncid,'geo_strf_velocity_mid_lat','double',[dims_test_cast_ml dims_test_cast_nl]);
netcdf.putAtt(ncid,geo_strf_velocity_mid_lat_ID,'computation_accuracy',gsw_cv.geo_strf_velocity_mid_lat_ca);

geo_strf_velocity_mid_long_ID = netcdf.defVar(ncid,'geo_strf_velocity_mid_long','double',[dims_test_cast_ml dims_test_cast_nl]);
netcdf.putAtt(ncid,geo_strf_velocity_mid_long_ID,'computation_accuracy',gsw_cv.geo_strf_velocity_mid_long_ca);

%% neutral versus isopycnal slopes and ratios

isopycnal_slope_ratio_ID = netcdf.defVar(ncid,'isopycnal_slope_ratio','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,isopycnal_slope_ratio_ID,'computation_accuracy',gsw_cv.isopycnal_slope_ratio_ca);

G_CT_ID = netcdf.defVar(ncid,'G_CT','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,G_CT_ID,'computation_accuracy',gsw_cv.G_CT_ca);

p_mid_G_CT_ID = netcdf.defVar(ncid,'p_mid_G_CT','double',[dims_test_cast_mm dims_test_cast_n]);
netcdf.putAtt(ncid,p_mid_G_CT_ID,'computation_accuracy',gsw_cv.p_mid_G_CT_ca);

ntpptCT_ID = netcdf.defVar(ncid,'ntpptCT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,ntpptCT_ID,'computation_accuracy',gsw_cv.ntpptCT_ca);

%% derivatives of enthalpy, entropy, CT and pt

CT_SA_ID = netcdf.defVar(ncid,'CT_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_SA_ID,'computation_accuracy',gsw_cv.CT_SA_ca);

CT_pt_ID = netcdf.defVar(ncid,'CT_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_pt_ID,'computation_accuracy',gsw_cv.CT_pt_ca);

CT_SA_SA_ID = netcdf.defVar(ncid,'CT_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_SA_SA_ID,'computation_accuracy',gsw_cv.CT_SA_SA_ca);

CT_SA_pt_ID = netcdf.defVar(ncid,'CT_SA_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_SA_pt_ID,'computation_accuracy',gsw_cv.CT_SA_pt_ca);

CT_pt_pt_ID = netcdf.defVar(ncid,'CT_pt_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_pt_pt_ID,'computation_accuracy',gsw_cv.CT_pt_pt_ca);

eta_SA_ID = netcdf.defVar(ncid,'eta_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,eta_SA_ID,'computation_accuracy',gsw_cv.eta_SA_ca);

eta_CT_ID = netcdf.defVar(ncid,'eta_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,eta_CT_ID,'computation_accuracy',gsw_cv.eta_CT_ca);

eta_SA_SA_ID = netcdf.defVar(ncid,'eta_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,eta_SA_SA_ID,'computation_accuracy',gsw_cv.eta_SA_SA_ca);

eta_SA_CT_ID = netcdf.defVar(ncid,'eta_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,eta_SA_CT_ID,'computation_accuracy',gsw_cv.eta_SA_CT_ca);

eta_CT_CT_ID = netcdf.defVar(ncid,'eta_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,eta_CT_CT_ID,'computation_accuracy',gsw_cv.eta_CT_CT_ca);

pt_SA_ID = netcdf.defVar(ncid,'pt_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_SA_ID,'computation_accuracy',gsw_cv.pt_SA_ca);

pt_CT_ID = netcdf.defVar(ncid,'pt_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_CT_ID,'computation_accuracy',gsw_cv.pt_CT_ca);

pt_SA_SA_ID = netcdf.defVar(ncid,'pt_SA_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_SA_SA_ID,'computation_accuracy',gsw_cv.pt_SA_SA_ca);

pt_SA_CT_ID = netcdf.defVar(ncid,'pt_SA_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_SA_CT_ID,'computation_accuracy',gsw_cv.pt_SA_CT_ca);

pt_CT_CT_ID = netcdf.defVar(ncid,'pt_CT_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pt_CT_CT_ID,'computation_accuracy',gsw_cv.pt_CT_CT_ca);

%% seawater properties at freezing temperatures

CT_freezing_ID = netcdf.defVar(ncid,'CT_freezing','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_freezing_ID,'computation_accuracy',gsw_cv.CT_freezing_ca);

CT_freezing_poly_ID = netcdf.defVar(ncid,'CT_freezing_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_freezing_poly_ID,'computation_accuracy',gsw_cv.CT_freezing_poly_ca);

t_freezing_ID = netcdf.defVar(ncid,'t_freezing','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_freezing_ID,'computation_accuracy',gsw_cv.t_freezing_ca);

t_freezing_poly_ID = netcdf.defVar(ncid,'t_freezing_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_freezing_poly_ID,'computation_accuracy',gsw_cv.t_freezing_poly_ca);

pot_enthalpy_ice_freezing_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_ca);

pot_enthalpy_ice_freezing_poly_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_poly_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_poly_ca);

SA_freezing_from_CT_ID = netcdf.defVar(ncid,'SA_freezing_from_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_freezing_from_CT_ID,'computation_accuracy',gsw_cv.SA_freezing_from_CT_ca);

SA_freezing_from_CT_poly_ID = netcdf.defVar(ncid,'SA_freezing_from_CT_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_freezing_from_CT_poly_ID,'computation_accuracy',gsw_cv.SA_freezing_from_CT_poly_ca);

SA_freezing_from_t_ID = netcdf.defVar(ncid,'SA_freezing_from_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_freezing_from_t_ID,'computation_accuracy',gsw_cv.SA_freezing_from_t_ca);

SA_freezing_from_t_poly_ID = netcdf.defVar(ncid,'SA_freezing_from_t_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_freezing_from_t_poly_ID,'computation_accuracy',gsw_cv.SA_freezing_from_t_poly_ca);

pressure_freezing_CT_ID = netcdf.defVar(ncid,'pressure_freezing_CT','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pressure_freezing_CT_ID,'computation_accuracy',gsw_cv.pressure_freezing_CT_ca);

CTfreezing_SA_ID = netcdf.defVar(ncid,'CTfreezing_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CTfreezing_SA_ID,'computation_accuracy',gsw_cv.CTfreezing_SA_ca);

CTfreezing_P_ID = netcdf.defVar(ncid,'CTfreezing_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CTfreezing_P_ID,'computation_accuracy',gsw_cv.CTfreezing_P_ca);

CTfreezing_SA_poly_ID = netcdf.defVar(ncid,'CTfreezing_SA_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CTfreezing_SA_poly_ID,'computation_accuracy',gsw_cv.CTfreezing_SA_poly_ca);

CTfreezing_P_poly_ID = netcdf.defVar(ncid,'CTfreezing_P_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CTfreezing_P_poly_ID,'computation_accuracy',gsw_cv.CTfreezing_P_poly_ca);

tfreezing_SA_ID = netcdf.defVar(ncid,'tfreezing_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,tfreezing_SA_ID,'computation_accuracy',gsw_cv.tfreezing_SA_ca);

tfreezing_P_ID = netcdf.defVar(ncid,'tfreezing_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,tfreezing_P_ID,'computation_accuracy',gsw_cv.tfreezing_P_ca);

tfreezing_SA_poly_ID = netcdf.defVar(ncid,'tfreezing_SA_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,tfreezing_SA_poly_ID,'computation_accuracy',gsw_cv.tfreezing_SA_poly_ca);

tfreezing_P_poly_ID = netcdf.defVar(ncid,'tfreezing_P_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,tfreezing_P_poly_ID,'computation_accuracy',gsw_cv.tfreezing_P_poly_ca);

pot_enthalpy_ice_freezing_SA_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing_SA','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_SA_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_SA_ca);

pot_enthalpy_ice_freezing_P_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing_P','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_P_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_P_ca);

pot_enthalpy_ice_freezing_SA_poly_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing_SA_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_SA_poly_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_SA_poly_ca);

pot_enthalpy_ice_freezing_P_poly_ID = netcdf.defVar(ncid,'pot_enthalpy_ice_freezing_P_poly','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_ice_freezing_P_poly_ID,'computation_accuracy',gsw_cv.pot_enthalpy_ice_freezing_P_poly_ca);

latentheat_melting_ID = netcdf.defVar(ncid,'latentheat_melting','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,latentheat_melting_ID,'computation_accuracy',gsw_cv.latentheat_melting_ca);

%%  thermodynamic interaction between ice and seawater
melting_ice_SA_CT_ratio_ID = netcdf.defVar(ncid,'melting_ice_SA_CT_ratio','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_SA_CT_ratio_ID,'computation_accuracy',gsw_cv.melting_ice_SA_CT_ratio_ca);

melting_ice_SA_CT_ratio_poly_ID = netcdf.defVar(ncid,'melting_ice_SA_CT_ratio_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_SA_CT_ratio_poly_ID,'computation_accuracy',gsw_cv.melting_ice_SA_CT_ratio_poly_ca);

melting_ice_equilibrium_SA_CT_ratio_ID = netcdf.defVar(ncid,'melting_ice_equilibrium_SA_CT_ratio','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_equilibrium_SA_CT_ratio_ID,'computation_accuracy',gsw_cv.melting_ice_equilibrium_SA_CT_ratio_ca);

melting_ice_equilibrium_SA_CT_ratio_poly_ID = netcdf.defVar(ncid,'melting_ice_equilibrium_SA_CT_ratio_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_equilibrium_SA_CT_ratio_poly_ID,'computation_accuracy',gsw_cv.melting_ice_equilibrium_SA_CT_ratio_poly_ca);

melting_ice_into_seawater_SA_final_ID = netcdf.defVar(ncid,'melting_ice_into_seawater_SA_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_into_seawater_SA_final_ID,'computation_accuracy',gsw_cv.melting_ice_into_seawater_SA_final_ca);

melting_ice_into_seawater_CT_final_ID = netcdf.defVar(ncid,'melting_ice_into_seawater_CT_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_ice_into_seawater_CT_final_ID,'computation_accuracy',gsw_cv.melting_ice_into_seawater_CT_final_ca);

ice_fraction_to_freeze_seawater_SA_freeze_ID = netcdf.defVar(ncid,'ice_fraction_to_freeze_seawater_SA_freeze','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,ice_fraction_to_freeze_seawater_SA_freeze_ID,'computation_accuracy',gsw_cv.ice_fraction_to_freeze_seawater_SA_freeze_ca);

ice_fraction_to_freeze_seawater_CT_freeze_ID = netcdf.defVar(ncid,'ice_fraction_to_freeze_seawater_CT_freeze','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,ice_fraction_to_freeze_seawater_CT_freeze_ID,'computation_accuracy',gsw_cv.ice_fraction_to_freeze_seawater_CT_freeze_ca);

ice_fraction_to_freeze_seawater_w_Ih_ID = netcdf.defVar(ncid,'ice_fraction_to_freeze_seawater_w_Ih','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,ice_fraction_to_freeze_seawater_w_Ih_ID,'computation_accuracy',gsw_cv.ice_fraction_to_freeze_seawater_w_Ih_ca);

dSA_dCT_frazil_ID = netcdf.defVar(ncid,'dSA_dCT_frazil','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dSA_dCT_frazil_ID,'computation_accuracy',gsw_cv.dSA_dCT_frazil_ca);

dSA_dP_frazil_ID = netcdf.defVar(ncid,'dSA_dP_frazil','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dSA_dP_frazil_ID,'computation_accuracy',gsw_cv.dSA_dP_frazil_ca);

dCT_dP_frazil_ID = netcdf.defVar(ncid,'dCT_dP_frazil','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dCT_dP_frazil_ID,'computation_accuracy',gsw_cv.dCT_dP_frazil_ca);

dSA_dCT_frazil_poly_ID = netcdf.defVar(ncid,'dSA_dCT_frazil_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dSA_dCT_frazil_poly_ID,'computation_accuracy',gsw_cv.dSA_dCT_frazil_poly_ca);

dSA_dP_frazil_poly_ID = netcdf.defVar(ncid,'dSA_dP_frazil_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dSA_dP_frazil_poly_ID,'computation_accuracy',gsw_cv.dSA_dP_frazil_poly_ca);

dCT_dP_frazil_poly_ID = netcdf.defVar(ncid,'dCT_dP_frazil_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,dCT_dP_frazil_poly_ID,'computation_accuracy',gsw_cv.dCT_dP_frazil_poly_ca);

frazil_properties_potential_SA_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_SA_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_SA_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_SA_final_ca);

frazil_properties_potential_CT_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_CT_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_CT_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_CT_final_ca);

frazil_properties_potential_w_Ih_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_w_Ih_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_w_Ih_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_w_Ih_final_ca);

frazil_properties_potential_poly_SA_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_poly_SA_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_poly_SA_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_poly_SA_final_ca);

frazil_properties_potential_poly_CT_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_poly_CT_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_poly_CT_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_poly_CT_final_ca);

frazil_properties_potential_poly_w_Ih_final_ID = netcdf.defVar(ncid,'frazil_properties_potential_poly_w_Ih_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_potential_poly_w_Ih_final_ID,'computation_accuracy',gsw_cv.frazil_properties_potential_poly_w_Ih_final_ca);

frazil_properties_SA_final_ID = netcdf.defVar(ncid,'frazil_properties_SA_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_SA_final_ID,'computation_accuracy',gsw_cv.frazil_properties_SA_final_ca);

frazil_properties_CT_final_ID = netcdf.defVar(ncid,'frazil_properties_CT_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_CT_final_ID,'computation_accuracy',gsw_cv.frazil_properties_CT_final_ca);

frazil_properties_w_Ih_final_ID = netcdf.defVar(ncid,'frazil_properties_w_Ih_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,frazil_properties_w_Ih_final_ID,'computation_accuracy',gsw_cv.frazil_properties_w_Ih_final_ca);

%%  thermodynamic interaction between sea ice and seawater

melting_seaice_SA_CT_ratio_ID = netcdf.defVar(ncid,'melting_seaice_SA_CT_ratio','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_SA_CT_ratio_ID,'computation_accuracy',gsw_cv.melting_seaice_SA_CT_ratio_ca);

melting_seaice_SA_CT_ratio_poly_ID = netcdf.defVar(ncid,'melting_seaice_SA_CT_ratio_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_SA_CT_ratio_poly_ID,'computation_accuracy',gsw_cv.melting_seaice_SA_CT_ratio_poly_ca);

melting_seaice_equilibrium_SA_CT_ratio_ID = netcdf.defVar(ncid,'melting_seaice_equilibrium_SA_CT_ratio','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_equilibrium_SA_CT_ratio_ID,'computation_accuracy',gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_ca);

melting_seaice_equilibrium_SA_CT_ratio_poly_ID = netcdf.defVar(ncid,'melting_seaice_equilibrium_SA_CT_ratio_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_equilibrium_SA_CT_ratio_poly_ID,'computation_accuracy',gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_poly_ca);

melting_seaice_into_seawater_SA_final_ID = netcdf.defVar(ncid,'melting_seaice_into_seawater_SA_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_into_seawater_SA_final_ID,'computation_accuracy',gsw_cv.melting_seaice_into_seawater_SA_final_ca);

melting_seaice_into_seawater_CT_final_ID = netcdf.defVar(ncid,'melting_seaice_into_seawater_CT_final','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,melting_seaice_into_seawater_CT_final_ID,'computation_accuracy',gsw_cv.melting_seaice_into_seawater_CT_final_ca);

seaice_fraction_to_freeze_seawater_SA_freeze_ID = netcdf.defVar(ncid,'seaice_fraction_to_freeze_seawater_SA_freeze','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,seaice_fraction_to_freeze_seawater_SA_freeze_ID,'computation_accuracy',gsw_cv.seaice_fraction_to_freeze_seawater_SA_freeze_ca);

seaice_fraction_to_freeze_seawater_CT_freeze_ID = netcdf.defVar(ncid,'seaice_fraction_to_freeze_seawater_CT_freeze','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,seaice_fraction_to_freeze_seawater_CT_freeze_ID,'computation_accuracy',gsw_cv.seaice_fraction_to_freeze_seawater_CT_freeze_ca);

seaice_fraction_to_freeze_seawater_w_Ih_ID = netcdf.defVar(ncid,'seaice_fraction_to_freeze_seawater_w_Ih','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,seaice_fraction_to_freeze_seawater_w_Ih_ID,'computation_accuracy',gsw_cv.seaice_fraction_to_freeze_seawater_w_Ih_ca);

%% themodynamic properties of ice Ih

specvol_ice_ID = netcdf.defVar(ncid,'specvol_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,specvol_ice_ID,'computation_accuracy',gsw_cv.specvol_ice_ca);

alpha_wrt_t_ice_ID = netcdf.defVar(ncid,'alpha_wrt_t_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,alpha_wrt_t_ice_ID,'computation_accuracy',gsw_cv.alpha_wrt_t_ice_ca);

rho_ice_ID = netcdf.defVar(ncid,'rho_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,rho_ice_ID,'computation_accuracy',gsw_cv.rho_ice_ca);

pressure_coefficient_ice_ID = netcdf.defVar(ncid,'pressure_coefficient_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pressure_coefficient_ice_ID,'computation_accuracy',gsw_cv.pressure_coefficient_ice_ca);

sound_speed_ice_ID = netcdf.defVar(ncid,'sound_speed_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,sound_speed_ice_ID,'computation_accuracy',gsw_cv.sound_speed_ice_ca);

kappa_ice_ID = netcdf.defVar(ncid,'kappa_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,kappa_ice_ID,'computation_accuracy',gsw_cv.kappa_ice_ca);

kappa_const_t_ice_ID = netcdf.defVar(ncid,'kappa_const_t_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,kappa_const_t_ice_ID,'computation_accuracy',gsw_cv.kappa_const_t_ice_ca);

internal_energy_ice_ID = netcdf.defVar(ncid,'internal_energy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,internal_energy_ice_ID,'computation_accuracy',gsw_cv.internal_energy_ice_ca);

enthalpy_ice_ID = netcdf.defVar(ncid,'enthalpy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_ice_ID,'computation_accuracy',gsw_cv.enthalpy_ice_ca);

entropy_ice_ID = netcdf.defVar(ncid,'entropy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,entropy_ice_ID,'computation_accuracy',gsw_cv.entropy_ice_ca);

cp_ice_ID = netcdf.defVar(ncid,'cp_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,cp_ice_ID,'computation_accuracy',gsw_cv.cp_ice_ca);

chem_potential_water_ice_ID = netcdf.defVar(ncid,'chem_potential_water_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,chem_potential_water_ice_ID,'computation_accuracy',gsw_cv.chem_potential_water_ice_ca);

Helmholtz_energy_ice_ID = netcdf.defVar(ncid,'Helmholtz_energy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,Helmholtz_energy_ice_ID,'computation_accuracy',gsw_cv.Helmholtz_energy_ice_ca);

adiabatic_lapse_rate_ice_ID = netcdf.defVar(ncid,'adiabatic_lapse_rate_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,adiabatic_lapse_rate_ice_ID,'computation_accuracy',gsw_cv.adiabatic_lapse_rate_ice_ca);

pt0_from_t_ice_ID = netcdf.defVar(ncid,'pt0_from_t_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pt0_from_t_ice_ID,'computation_accuracy',gsw_cv.pt0_from_t_ice_ca);

pt_from_t_ice_ID = netcdf.defVar(ncid,'pt_from_t_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pt_from_t_ice_ID,'computation_accuracy',gsw_cv.pt_from_t_ice_ca);

t_from_pt0_ice_ID = netcdf.defVar(ncid,'t_from_pt0_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,t_from_pt0_ice_ID,'computation_accuracy',gsw_cv.t_from_pt0_ice_ca);

t_from_rho_ice_ID = netcdf.defVar(ncid,'t_from_rho_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,t_from_rho_ice_ID,'computation_accuracy',gsw_cv.t_from_rho_ice_ca);

pot_enthalpy_from_pt_ice_ID = netcdf.defVar(ncid,'pot_enthalpy_from_pt_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_from_pt_ice_ID,'computation_accuracy',gsw_cv.pot_enthalpy_from_pt_ice_ca);

pt_from_pot_enthalpy_ice_ID = netcdf.defVar(ncid,'pt_from_pot_enthalpy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pt_from_pot_enthalpy_ice_ID,'computation_accuracy',gsw_cv.pt_from_pot_enthalpy_ice_ca);

pot_enthalpy_from_pt_ice_poly_ID = netcdf.defVar(ncid,'pot_enthalpy_from_pt_ice_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_from_pt_ice_poly_ID,'computation_accuracy',gsw_cv.pot_enthalpy_from_pt_ice_poly_ca);

pt_from_pot_enthalpy_ice_poly_ID = netcdf.defVar(ncid,'pt_from_pot_enthalpy_ice_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pt_from_pot_enthalpy_ice_poly_ID,'computation_accuracy',gsw_cv.pt_from_pot_enthalpy_ice_poly_ca);

pot_enthalpy_from_specvol_ice_ID = netcdf.defVar(ncid,'pot_enthalpy_from_specvol_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_from_specvol_ice_ID,'computation_accuracy',gsw_cv.pot_enthalpy_from_specvol_ice_ca);

specvol_from_pot_enthalpy_ice_ID = netcdf.defVar(ncid,'specvol_from_pot_enthalpy_ice','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,specvol_from_pot_enthalpy_ice_ID,'computation_accuracy',gsw_cv.specvol_from_pot_enthalpy_ice_ca);

pot_enthalpy_from_specvol_ice_poly_ID = netcdf.defVar(ncid,'pot_enthalpy_from_specvol_ice_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,pot_enthalpy_from_specvol_ice_poly_ID,'computation_accuracy',gsw_cv.pot_enthalpy_from_specvol_ice_poly_ca);

specvol_from_pot_enthalpy_ice_poly_ID = netcdf.defVar(ncid,'specvol_from_pot_enthalpy_ice_poly','double',[dims_Arctic_test_cast_m dims_Arctic_test_cast_n]);
netcdf.putAtt(ncid,specvol_from_pot_enthalpy_ice_poly_ID,'computation_accuracy',gsw_cv.specvol_from_pot_enthalpy_ice_poly_ca);

%% isobaric evaporation enthalpy 

latentheat_evap_CT_ID = netcdf.defVar(ncid,'latentheat_evap_CT','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,latentheat_evap_CT_ID,'computation_accuracy',gsw_cv.latentheat_evap_CT_ca);

latentheat_evap_t_ID = netcdf.defVar(ncid,'latentheat_evap_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,latentheat_evap_t_ID,'computation_accuracy',gsw_cv.latentheat_evap_t_ca);

%% spiciness

spiciness0_ID = netcdf.defVar(ncid,'spiciness0','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,spiciness0_ID,'computation_accuracy',gsw_cv.spiciness0_ca);

spiciness1_ID = netcdf.defVar(ncid,'spiciness1','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,spiciness1_ID,'computation_accuracy',gsw_cv.spiciness1_ca);

spiciness2_ID = netcdf.defVar(ncid,'spiciness2','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,spiciness2_ID,'computation_accuracy',gsw_cv.spiciness2_ca);

SA_from_sigma0_spiciness0_ID = netcdf.defVar(ncid,'SA_from_sigma0_spiciness0','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_sigma0_spiciness0_ID,'computation_accuracy',gsw_cv.SA_spicsig0_ca);

CT_from_sigma0_spiciness0_ID = netcdf.defVar(ncid,'CT_from_sigma0_spiciness0','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_sigma0_spiciness0_ID,'computation_accuracy',gsw_cv.CT_spicsig0_ca);

SA_from_sigma1_spiciness1_ID = netcdf.defVar(ncid,'SA_from_sigma1_spiciness1','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_sigma1_spiciness1_ID,'computation_accuracy',gsw_cv.SA_spicsig1_ca);

CT_from_sigma1_spiciness1_ID = netcdf.defVar(ncid,'CT_from_sigma1_spiciness1','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_sigma1_spiciness1_ID,'computation_accuracy',gsw_cv.CT_spicsig1_ca);

SA_from_sigma2_spiciness2_ID = netcdf.defVar(ncid,'SA_from_sigma2_spiciness2','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_sigma2_spiciness2_ID,'computation_accuracy',gsw_cv.SA_spicsig2_ca);

CT_from_sigma2_spiciness2_ID = netcdf.defVar(ncid,'CT_from_sigma2_spiciness2','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_sigma2_spiciness2_ID,'computation_accuracy',gsw_cv.CT_spicsig2_ca);


%% planet earth properties

f_ID = netcdf.defVar(ncid,'f','double',dims_test_cast_n);
netcdf.putAtt(ncid,f_ID,'computation_accuracy',gsw_cv.f_ca);

grav_ID = netcdf.defVar(ncid,'grav','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,grav_ID,'computation_accuracy',gsw_cv.grav_ca);

distance_ID = netcdf.defVar(ncid,'distance','double', [dims_test_cast_ml dims_test_cast_nl]);
netcdf.putAtt(ncid,distance_ID,'computation_accuracy',gsw_cv.distance_ca);

%% TEOS-10 constants
T0_ID = netcdf.defVar(ncid,'T0','double',dims_value_test_cast);

P0_ID = netcdf.defVar(ncid,'P0','double', dims_value_test_cast);

SSO_ID = netcdf.defVar(ncid,'SSO','double', dims_value_test_cast);

uPS_ID = netcdf.defVar(ncid,'uPS','double',dims_value_test_cast);

cp0_ID = netcdf.defVar(ncid,'cp0','double', dims_value_test_cast);

C3515_ID = netcdf.defVar(ncid,'C3515','double', dims_value_test_cast);

SonCl_ID = netcdf.defVar(ncid,'SonCl','double', dims_value_test_cast);

valence_factor_ID = netcdf.defVar(ncid,'valence_factor','double', dims_value_test_cast);

atomic_weight_ID = netcdf.defVar(ncid,'atomic_weight','double', dims_value_test_cast);

%% dissolved gasses

Arsol_ID = netcdf.defVar(ncid,'Arsol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Arsol_ID,'computation_accuracy',gsw_cv.Arsol_ca);

Arsol_SP_pt_ID = netcdf.defVar(ncid,'Arsol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Arsol_SP_pt_ID,'computation_accuracy',gsw_cv.Arsol_SP_pt_ca);

Hesol_ID = netcdf.defVar(ncid,'Hesol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Hesol_ID,'computation_accuracy',gsw_cv.Hesol_ca);

Hesol_SP_pt_ID = netcdf.defVar(ncid,'Hesol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Hesol_SP_pt_ID,'computation_accuracy',gsw_cv.Hesol_SP_pt_ca);

Krsol_ID = netcdf.defVar(ncid,'Krsol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Krsol_ID,'computation_accuracy',gsw_cv.Krsol_ca);

Krsol_SP_pt_ID = netcdf.defVar(ncid,'Krsol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Krsol_SP_pt_ID,'computation_accuracy',gsw_cv.Krsol_SP_pt_ca);

N2Osol_ID = netcdf.defVar(ncid,'N2Osol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,N2Osol_ID,'computation_accuracy',gsw_cv.N2Osol_ca);

N2Osol_SP_pt_ID = netcdf.defVar(ncid,'N2Osol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,N2Osol_SP_pt_ID,'computation_accuracy',gsw_cv.N2Osol_SP_pt_ca);

N2sol_ID = netcdf.defVar(ncid,'N2sol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,N2sol_ID,'computation_accuracy',gsw_cv.N2sol_ca);

N2sol_SP_pt_ID = netcdf.defVar(ncid,'N2sol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,N2sol_SP_pt_ID,'computation_accuracy',gsw_cv.N2sol_SP_pt_ca);

Nesol_ID = netcdf.defVar(ncid,'Nesol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Nesol_ID,'computation_accuracy',gsw_cv.Nesol_ca);

Nesol_SP_pt_ID = netcdf.defVar(ncid,'Nesol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Nesol_SP_pt_ID,'computation_accuracy',gsw_cv.Nesol_SP_pt_ca);

O2sol_ID = netcdf.defVar(ncid,'O2sol','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,O2sol_ID,'computation_accuracy',gsw_cv.O2sol_ca);

O2sol_SP_pt_ID = netcdf.defVar(ncid,'O2sol_SP_pt','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,O2sol_SP_pt_ID,'computation_accuracy',gsw_cv.O2sol_SP_pt_ca);

%% density and enthalpy in terms of CT, derived from the exact Gibbs function

specvol_CT_exact_ID = netcdf.defVar(ncid,'specvol_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_CT_exact_ID,'computation_accuracy',gsw_cv.specvol_CT_exact_ca);

alpha_CT_exact_ID = netcdf.defVar(ncid,'alpha_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_CT_exact_ID,'computation_accuracy',gsw_cv.alpha_CT_exact_ca);

beta_CT_exact_ID = netcdf.defVar(ncid,'beta_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_CT_exact_ID,'computation_accuracy',gsw_cv.beta_CT_exact_ca);

alpha_on_beta_CT_exact_ID = netcdf.defVar(ncid,'alpha_on_beta_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_on_beta_CT_exact_ID,'computation_accuracy',gsw_cv.alpha_on_beta_CT_exact_ca);

v_vab_CT_exact_ID = netcdf.defVar(ncid,'v_vab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_vab_CT_exact_ID,'computation_accuracy',gsw_cv.v_vab_CT_exact_ca);

alpha_vab_CT_exact_ID = netcdf.defVar(ncid,'alpha_vab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_vab_CT_exact_ID,'computation_accuracy',gsw_cv.alpha_vab_CT_exact_ca);

beta_vab_CT_exact_ID = netcdf.defVar(ncid,'beta_vab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_vab_CT_exact_ID,'computation_accuracy',gsw_cv.beta_vab_CT_exact_ca);

v_SA_CT_exact_ID = netcdf.defVar(ncid,'v_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_CT_exact_ca);

v_CT_CT_exact_ID = netcdf.defVar(ncid,'v_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_CT_exact_ID,'computation_accuracy',gsw_cv.v_CT_CT_exact_ca);

v_P_CT_exact_ID = netcdf.defVar(ncid,'v_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_P_CT_exact_ID,'computation_accuracy',gsw_cv.v_P_CT_exact_ca);

v_SA_SA_CT_exact_ID = netcdf.defVar(ncid,'v_SA_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_SA_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_SA_CT_exact_ca);

v_SA_CT_CT_exact_ID = netcdf.defVar(ncid,'v_SA_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_CT_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_CT_CT_exact_ca);

v_CT_CT_CT_exact_ID = netcdf.defVar(ncid,'v_CT_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_CT_CT_exact_ID,'computation_accuracy',gsw_cv.v_CT_CT_CT_exact_ca);

v_SA_P_CT_exact_ID = netcdf.defVar(ncid,'v_SA_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_P_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_P_CT_exact_ca);

v_CT_P_CT_exact_ID = netcdf.defVar(ncid,'v_CT_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_CT_P_CT_exact_ID,'computation_accuracy',gsw_cv.v_CT_P_CT_exact_ca);

v_SA_wrt_h_CT_exact_ID = netcdf.defVar(ncid,'v_SA_wrt_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_wrt_h_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_wrt_h_CT_exact_ca);

v_h_CT_exact_ID = netcdf.defVar(ncid,'v_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_h_CT_exact_ID,'computation_accuracy',gsw_cv.v_h_CT_exact_ca);

v_SA_SA_wrt_h_CT_exact_ID = netcdf.defVar(ncid,'v_SA_SA_wrt_h_CT_exact_ID','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_SA_wrt_h_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_SA_wrt_h_CT_exact_ca);

v_SA_h_CT_exact_ID = netcdf.defVar(ncid,'v_SA_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_SA_h_CT_exact_ID,'computation_accuracy',gsw_cv.v_SA_h_CT_exact_ca);

v_h_h_CT_exact_ID = netcdf.defVar(ncid,'v_h_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,v_h_h_CT_exact_ID,'computation_accuracy',gsw_cv.v_h_h_CT_exact_ca);

specvol_anom_CT_exact_ID = netcdf.defVar(ncid,'specvol_anom_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_anom_CT_exact_ID,'computation_accuracy',gsw_cv.specvol_anom_CT_exact_ca);

specvol_anom_standard_CT_exact_ID = netcdf.defVar(ncid,'specvol_anom_standard_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_anom_CT_exact_ID,'computation_accuracy',gsw_cv.specvol_anom_standard_CT_exact_ca);

alpha_rab_CT_exact_ID = netcdf.defVar(ncid,'alpha_rab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);

beta_rab_CT_exact_ID = netcdf.defVar(ncid,'beta_rab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);

rho_rab_CT_exact_ID = netcdf.defVar(ncid,'rho_rab_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);

rho_CT_exact_ID = netcdf.defVar(ncid,'rho_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_exact_ID,'computation_accuracy',gsw_cv.rho_CT_exact_ca);

rho_SA_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_CT_exact_ca);

rho_CT_CT_exact_ID = netcdf.defVar(ncid,'rho_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_CT_exact_ID,'computation_accuracy',gsw_cv.rho_CT_CT_exact_ca);

rho_P_CT_exact_ID = netcdf.defVar(ncid,'rho_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_P_CT_exact_ID,'computation_accuracy',gsw_cv.rho_P_CT_exact_ca);

rho_SA_SA_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_SA_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_SA_CT_exact_ca);

rho_SA_CT_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_CT_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_CT_CT_exact_ca);

rho_CT_CT_CT_exact_ID = netcdf.defVar(ncid,'rho_CT_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_CT_CT_exact_ID,'computation_accuracy',gsw_cv.rho_CT_CT_CT_exact_ca);

rho_SA_P_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_P_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_P_CT_exact_ca);

rho_CT_P_CT_exact_ID = netcdf.defVar(ncid,'rho_CT_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_CT_P_CT_exact_ID,'computation_accuracy',gsw_cv.rho_CT_P_CT_exact_ca);

rho_SA_wrt_h_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_wrt_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_wrt_h_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_wrt_h_CT_exact_ca);

rho_h_CT_exact_ID = netcdf.defVar(ncid,'rho_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_h_CT_exact_ID,'computation_accuracy',gsw_cv.rho_h_CT_exact_ca);

rho_SA_SA_wrt_h_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_SA_wrt_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_SA_wrt_h_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_SA_wrt_h_CT_exact_ca);

rho_SA_h_CT_exact_ID = netcdf.defVar(ncid,'rho_SA_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_SA_h_CT_exact_ID,'computation_accuracy',gsw_cv.rho_SA_h_CT_exact_ca);

rho_h_h_CT_exact_ID = netcdf.defVar(ncid,'rho_h_h_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_h_h_CT_exact_ID,'computation_accuracy',gsw_cv.rho_h_h_CT_exact_ca);

sigma0_CT_exact_ID = netcdf.defVar(ncid,'sigma0_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma0_CT_exact_ID,'computation_accuracy',gsw_cv.sigma0_CT_exact_ca);

sigma1_CT_exact_ID = netcdf.defVar(ncid,'sigma1_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma1_CT_exact_ID,'computation_accuracy',gsw_cv.sigma1_CT_exact_ca);

sigma2_CT_exact_ID = netcdf.defVar(ncid,'sigma2_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma2_CT_exact_ID,'computation_accuracy',gsw_cv.sigma2_CT_exact_ca);

sigma3_CT_exact_ID = netcdf.defVar(ncid,'sigma3_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma3_CT_exact_ID,'computation_accuracy',gsw_cv.sigma3_CT_exact_ca);

sigma4_CT_exact_ID = netcdf.defVar(ncid,'sigma4_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma4_CT_exact_ID,'computation_accuracy',gsw_cv.sigma4_CT_exact_ca);

cabbeling_CT_exact_ID = netcdf.defVar(ncid,'cabbeling_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,cabbeling_CT_exact_ID,'computation_accuracy',gsw_cv.cabbeling_CT_exact_ca);

thermobaric_CT_exact_ID = netcdf.defVar(ncid,'thermobaric_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,thermobaric_CT_exact_ID,'computation_accuracy',gsw_cv.thermobaric_CT_exact_ca);

enthalpy_CT_exact_ID = netcdf.defVar(ncid,'enthalpy_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_CT_exact_ID,'computation_accuracy',gsw_cv.enthalpy_CT_exact_ca);

enthalpy_diff_CT_exact_ID = netcdf.defVar(ncid,'enthalpy_diff_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_diff_CT_exact_ID,'computation_accuracy',gsw_cv.enthalpy_diff_CT_exact_ca);

dynamic_enthalpy_CT_exact_ID = netcdf.defVar(ncid,'dynamic_enthalpy_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,dynamic_enthalpy_CT_exact_ID,'computation_accuracy',gsw_cv.dynamic_enthalpy_CT_exact_ca);

h_SA_CT_exact_ID = netcdf.defVar(ncid,'h_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_CT_exact_ID,'computation_accuracy',gsw_cv.h_SA_CT_exact_ca);

h_CT_CT_exact_ID = netcdf.defVar(ncid,'h_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_CT_CT_exact_ID,'computation_accuracy',gsw_cv.h_CT_CT_exact_ca);

h_SA_SA_CT_exact_ID = netcdf.defVar(ncid,'h_SA_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_SA_CT_exact_ID,'computation_accuracy',gsw_cv.h_SA_SA_CT_exact_ca);

h_SA_CT_CT_exact_ID = netcdf.defVar(ncid,'h_SA_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_CT_CT_exact_ID,'computation_accuracy',gsw_cv.h_SA_CT_CT_exact_ca);

h_CT_CT_CT_exact_ID = netcdf.defVar(ncid,'h_CT_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_CT_CT_CT_exact_ID,'computation_accuracy',gsw_cv.h_CT_CT_CT_exact_ca);

sound_speed_CT_exact_ID = netcdf.defVar(ncid,'sound_speed_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sound_speed_CT_exact_ID,'computation_accuracy',gsw_cv.sound_speed_CT_exact_ca);

kappa_CT_exact_ID = netcdf.defVar(ncid,'kappa_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,kappa_CT_exact_ID,'computation_accuracy',gsw_cv.kappa_CT_exact_ca);

internal_energy_CT_exact_ID = netcdf.defVar(ncid,'internal_energy_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,internal_energy_CT_exact_ID,'computation_accuracy',gsw_cv.internal_energy_CT_exact_ca);

u_SA_CT_exact_ID = netcdf.defVar(ncid,'u_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_CT_exact_ID,'computation_accuracy',gsw_cv.u_SA_CT_exact_ca);

u_CT_CT_exact_ID = netcdf.defVar(ncid,'u_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_CT_exact_ID,'computation_accuracy',gsw_cv.u_CT_CT_exact_ca);

u_P_CT_exact_ID = netcdf.defVar(ncid,'u_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_P_CT_exact_ID,'computation_accuracy',gsw_cv.u_P_CT_exact_ca);

u_SA_SA_CT_exact_ID = netcdf.defVar(ncid,'u_SA_SA_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_SA_CT_exact_ID,'computation_accuracy',gsw_cv.u_SA_SA_CT_exact_ca);

u_SA_CT_CT_exact_ID = netcdf.defVar(ncid,'u_SA_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_CT_CT_exact_ID,'computation_accuracy',gsw_cv.u_SA_CT_CT_exact_ca);

u_CT_CT_CT_exact_ID = netcdf.defVar(ncid,'u_CT_CT_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_CT_CT_exact_ID,'computation_accuracy',gsw_cv.u_CT_CT_CT_exact_ca);

u_SA_P_CT_exact_ID = netcdf.defVar(ncid,'u_SA_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_SA_P_CT_exact_ID,'computation_accuracy',gsw_cv.u_SA_P_CT_exact_ca);

u_CT_P_CT_exact_ID = netcdf.defVar(ncid,'u_CT_P_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,u_CT_P_CT_exact_ID,'computation_accuracy',gsw_cv.u_CT_P_CT_exact_ca);

CT_from_enthalpy_exact_ID = netcdf.defVar(ncid,'CT_from_enthalpy_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_enthalpy_exact_ID,'computation_accuracy',gsw_cv.CT_from_enthalpy_exact_ca);

SA_from_rho_CT_exact_ID = netcdf.defVar(ncid,'SA_from_rho_CT_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_rho_CT_exact_ID,'computation_accuracy',gsw_cv.SA_from_rho_CT_exact_ca);

CT_from_rho_exact_ID = netcdf.defVar(ncid,'CT_from_rho_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_from_rho_exact_ID,'computation_accuracy',gsw_cv.CT_from_rho_exact_ca);

CT_maxdensity_exact_ID = netcdf.defVar(ncid,'CT_maxdensity_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_maxdensity_exact_ID,'computation_accuracy',gsw_cv.CT_maxdensity_exact_ca);

%% Labrortory functions

rho_t_exact_ID = netcdf.defVar(ncid,'rho_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,rho_t_exact_ID,'computation_accuracy',gsw_cv.rho_t_exact_ca);

SA_from_rho_t_exact_ID = netcdf.defVar(ncid,'SA_from_rho_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,SA_from_rho_t_exact_ID,'computation_accuracy',gsw_cv.SA_from_rho_t_exact_ca);

deltaSA_from_rho_t_exact_ID = netcdf.defVar(ncid,'deltaSA_from_rho_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,deltaSA_from_rho_t_exact_ID,'computation_accuracy',gsw_cv.deltaSA_from_rho_t_exact_ca);

%% basic thermodynamic properties interms of in-situ t, derived from the exact Gibbs function

specvol_t_exact_ID = netcdf.defVar(ncid,'specvol_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_t_exact_ID,'computation_accuracy',gsw_cv.specvol_t_exact_ca);

alpha_wrt_CT_t_exact_ID = netcdf.defVar(ncid,'alpha_wrt_CT_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_wrt_CT_t_exact_ID,'computation_accuracy',gsw_cv.alpha_wrt_CT_t_exact_ca);

alpha_wrt_pt_t_exact_ID = netcdf.defVar(ncid,'alpha_wrt_pt_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_wrt_pt_t_exact_ID,'computation_accuracy',gsw_cv.alpha_wrt_pt_t_exact_ca);

alpha_wrt_t_exact_ID = netcdf.defVar(ncid,'alpha_wrt_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,alpha_wrt_t_exact_ID,'computation_accuracy',gsw_cv.alpha_wrt_t_exact_ca);

beta_const_CT_t_exact_ID = netcdf.defVar(ncid,'beta_const_CT_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_const_CT_t_exact_ID,'computation_accuracy',gsw_cv.beta_const_CT_t_exact_ca);

beta_const_pt_t_exact_ID = netcdf.defVar(ncid,'beta_const_pt_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_const_pt_t_exact_ID,'computation_accuracy',gsw_cv.beta_const_pt_t_exact_ca);

beta_const_t_exact_ID = netcdf.defVar(ncid,'beta_const_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,beta_const_t_exact_ID,'computation_accuracy',gsw_cv.beta_const_t_exact_ca);

specvol_anom_standard_t_exact_ID = netcdf.defVar(ncid,'specvol_anom_standard_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,specvol_anom_standard_t_exact_ID,'computation_accuracy',gsw_cv.specvol_anom_standard_t_exact_ca);

pot_rho_t_exact_ID = netcdf.defVar(ncid,'pot_rho_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,pot_rho_t_exact_ID,'computation_accuracy',gsw_cv.pot_rho_t_exact_ca);

sigma0_pt0_exact_ID = netcdf.defVar(ncid,'sigma0_pt0_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sigma0_pt0_exact_ID,'computation_accuracy',gsw_cv.sigma0_pt0_exact_ca);

enthalpy_t_exact_ID = netcdf.defVar(ncid,'enthalpy_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,enthalpy_t_exact_ID,'computation_accuracy',gsw_cv.enthalpy_t_exact_ca);

dynamic_enthalpy_t_exact_ID = netcdf.defVar(ncid,'dynamic_enthalpy_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,dynamic_enthalpy_t_exact_ID,'computation_accuracy',gsw_cv.dynamic_enthalpy_t_exact_ca);

CT_SA_wrt_t_ID = netcdf.defVar(ncid,'CT_SA_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_SA_wrt_t_ID,'computation_accuracy',gsw_cv.CT_SA_wrt_t_ca);

CT_T_wrt_t_ID = netcdf.defVar(ncid,'CT_T_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_T_wrt_t_ID,'computation_accuracy',gsw_cv.CT_T_wrt_t_ca);

CT_P_wrt_t_ID = netcdf.defVar(ncid,'CT_P_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,CT_P_wrt_t_ID,'computation_accuracy',gsw_cv.CT_P_wrt_t_ca);

h_SA_wrt_t_ID = netcdf.defVar(ncid,'h_SA_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_SA_wrt_t_ID,'computation_accuracy',gsw_cv.h_SA_wrt_t_ca);

h_T_wrt_t_ID = netcdf.defVar(ncid,'h_T_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_T_wrt_t_ID,'computation_accuracy',gsw_cv.h_T_wrt_t_ca);

h_P_wrt_t_ID = netcdf.defVar(ncid,'h_P_wrt_t','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,h_P_wrt_t_ID,'computation_accuracy',gsw_cv.h_P_wrt_t_ca);

sound_speed_t_exact_ID = netcdf.defVar(ncid,'sound_speed_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,sound_speed_t_exact_ID,'computation_accuracy',gsw_cv.sound_speed_t_exact_ca);

kappa_t_exact_ID = netcdf.defVar(ncid,'kappa_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,kappa_t_exact_ID,'computation_accuracy',gsw_cv.kappa_t_exact_ca);

kappa_const_t_exact_ID = netcdf.defVar(ncid,'kappa_const_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,kappa_const_t_exact_ID,'computation_accuracy',gsw_cv.kappa_const_t_exact_ca);

internal_energy_t_exact_ID = netcdf.defVar(ncid,'internal_energy_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,internal_energy_t_exact_ID,'computation_accuracy',gsw_cv.internal_energy_t_exact_ca);

t_from_rho_exact_ID = netcdf.defVar(ncid,'t_from_rho_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_from_rho_exact_ID,'computation_accuracy',gsw_cv.t_from_rho_exact_ca);

t_maxdensity_exact_ID = netcdf.defVar(ncid,'t_maxdensity_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_maxdensity_exact_ID,'computation_accuracy',gsw_cv.t_maxdensity_exact_ca);

cp_t_exact_ID = netcdf.defVar(ncid,'cp_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,cp_t_exact_ID,'computation_accuracy',gsw_cv.cp_t_exact_ca);

isochoric_heat_cap_t_exact_ID = netcdf.defVar(ncid,'isochoric_heat_cap_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,isochoric_heat_cap_t_exact_ID,'computation_accuracy',gsw_cv.isochoric_heat_cap_t_exact_ca);

chem_potential_relative_t_exact_ID = netcdf.defVar(ncid,'chem_potential_relative_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,chem_potential_relative_t_exact_ID,'computation_accuracy',gsw_cv.chem_potential_relative_t_exact_ca);

chem_potential_water_t_exact_ID = netcdf.defVar(ncid,'chem_potential_water_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,chem_potential_water_t_exact_ID,'computation_accuracy',gsw_cv.chem_potential_water_t_exact_ca);

chem_potential_salt_t_exact_ID = netcdf.defVar(ncid,'chem_potential_salt_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,chem_potential_salt_t_exact_ID,'computation_accuracy',gsw_cv.chem_potential_salt_t_exact_ca);

t_deriv_chem_potential_water_t_exact_ID = netcdf.defVar(ncid,'t_deriv_chem_potential_water_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,t_deriv_chem_potential_water_t_exact_ID,'computation_accuracy',gsw_cv.t_deriv_chem_potential_water_t_exact_ca);

dilution_coefficient_t_exact_ID = netcdf.defVar(ncid,'dilution_coefficient_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,dilution_coefficient_t_exact_ID,'computation_accuracy',gsw_cv.dilution_coefficient_t_exact_ca);

Gibbs_energy_t_exact_ID = netcdf.defVar(ncid,'Gibbs_energy_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Gibbs_energy_t_exact_ID,'computation_accuracy',gsw_cv.Gibbs_energy_t_exact_ca);

Helmholtz_energy_t_exact_ID = netcdf.defVar(ncid,'Helmholtz_energy_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Helmholtz_energy_t_exact_ID,'computation_accuracy',gsw_cv.Helmholtz_energy_t_exact_ca);

osmotic_coefficient_t_exact_ID = netcdf.defVar(ncid,'osmotic_coefficient_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,osmotic_coefficient_t_exact_ID,'computation_accuracy',gsw_cv.osmotic_coefficient_t_exact_ca);

osmotic_pressure_t_exact_ID = netcdf.defVar(ncid,'osmotic_pressure_t_exact','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,osmotic_pressure_t_exact_ID,'computation_accuracy',gsw_cv.osmotic_pressure_t_exact_ca);

%% Library

Fdelta_ID = netcdf.defVar(ncid,'Fdelta','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,Fdelta_ID,'computation_accuracy',gsw_cv.Fdelta_ca);

deltaSA_atlas_ID = netcdf.defVar(ncid,'deltaSA_atlas','double',[dims_test_cast_m dims_test_cast_n]);
netcdf.putAtt(ncid,deltaSA_atlas_ID,'computation_accuracy',gsw_cv.deltaSA_atlas_ca);

%=============================================================================

%We are done defining the NetCdf
netcdf.endDef(ncid);
 
%Then store the dimension variables in
netcdf.putVar(ncid,p_ref_ID,p_ref);
netcdf.putVar(ncid,lats_ref_ID,lats_ref);
netcdf.putVar(ncid,longs_ref_ID,longs_ref);
netcdf.putVar(ncid,p_ref_cast_ID,p_ref_cast);
netcdf.putVar(ncid,lat_ref_cast_ID,lat_ref_cast);
netcdf.putVar(ncid,long_ref_cast_ID,long_ref_cast);

%Then store the main variables
netcdf.putVar(ncid,SAAR_ref_ID,SAAR_ref);
netcdf.putVar(ncid,SR_ref_ID,SR_ref);
netcdf.putVar(ncid,deltaSA_ref_ID,deltaSA_ref);
netcdf.putVar(ncid,ocean_ref_ID,ocean_ref);
netcdf.putVar(ncid,ndepth_ref_ID,ndepth_ref);
netcdf.putVar(ncid,SA_ref_cast_ID,SA_ref_cast);
netcdf.putVar(ncid,CT_ref_cast_ID,CT_ref_cast);
netcdf.putVar(ncid,gamma_n_ref_cast_ID,gamma_n_ref_cast);
netcdf.putVar(ncid,sigma_2_ref_cast_ID,sigma_2_ref_cast);

netcdf.putVar(ncid,SA_Arctic_ID,gsw_cv.SA_Arctic);
netcdf.putVar(ncid,t_Arctic_ID,gsw_cv.t_Arctic);
netcdf.putVar(ncid,p_Arctic_ID,gsw_cv.p_Arctic);
netcdf.putVar(ncid,latitude_Arctic_ID,gsw_cv.latitude_Arctic);
netcdf.putVar(ncid,longitude_Arctic_ID,gsw_cv.longitude_Arctic);
netcdf.putVar(ncid,SA_seaice_ID,gsw_cv.SA_seaice);
netcdf.putVar(ncid,t_seaice_ID,gsw_cv.t_seaice);
netcdf.putVar(ncid,w_seaice_ID,gsw_cv.w_seaice);
netcdf.putVar(ncid,CT_Arctic_ID,gsw_cv.CT_Arctic);
netcdf.putVar(ncid,t_ice_ID,gsw_cv.t_ice);
netcdf.putVar(ncid,w_ice_ID,gsw_cv.w_ice);
netcdf.putVar(ncid,SA_bulk_ID,gsw_cv.SA_bulk);
netcdf.putVar(ncid,h_bulk_ID,gsw_cv.h_bulk);
netcdf.putVar(ncid,h_pot_bulk_ID,gsw_cv.h_pot_bulk);
netcdf.putVar(ncid,SP_chck_cast_ID,gsw_cv.SP_chck_cast);
netcdf.putVar(ncid,t_chck_cast_ID,gsw_cv.t_chck_cast);
netcdf.putVar(ncid,p_chck_cast_ID,gsw_cv.p_chck_cast);
netcdf.putVar(ncid,lat_chck_cast_ID,gsw_cv.lat_chck_cast);
netcdf.putVar(ncid,long_chck_cast_ID,gsw_cv.long_chck_cast);
netcdf.putVar(ncid,pr_ID,gsw_cv.pr);
netcdf.putVar(ncid,pr_05_ID,gsw_cv.pr_05);
netcdf.putVar(ncid,p_chck_cast_shallow_ID,gsw_cv.p_chck_cast_shallow);
netcdf.putVar(ncid,p_chck_cast_deep_ID,gsw_cv.p_chck_cast_deep);
netcdf.putVar(ncid,delta_p_chck_cast_ID,gsw_cv.delta_p_chck_cast);
netcdf.putVar(ncid,Neutral_Density_ID,gsw_cv.Neutral_Density);
netcdf.putVar(ncid,p_Neutral_Density_ID,gsw_cv.p_Neutral_Density);

%% Practical Salinity (SP):- PSS-78

netcdf.putVar(ncid,C_from_SP_ID,gsw_cv.C_from_SP);
netcdf.putVar(ncid,SP_from_C_ID,gsw_cv.SP_from_C);

netcdf.putVar(ncid,R_from_SP_ID,gsw_cv.R_from_SP);
netcdf.putVar(ncid,SP_from_R_ID,gsw_cv.SP_from_R);

netcdf.putVar(ncid,Rt_chck_cast_ID,gsw_cv.Rt_chck_cast);

netcdf.putVar(ncid,SP_salinometer_ID,gsw_cv.SP_salinometer);

netcdf.putVar(ncid,SK_chck_cast_ID,gsw_cv.SK_chck_cast);
netcdf.putVar(ncid,SP_from_SK_ID,gsw_cv.SP_from_SK);

%% Absolute Salinity (SA), Preformed Salinity (Sstar) and Conservative Temperature (CT)

netcdf.putVar(ncid,SA_chck_cast_ID,gsw_cv.SA_chck_cast);
netcdf.putVar(ncid,SA_from_SP_ID,gsw_cv.SA_from_SP);
netcdf.putVar(ncid,Sstar_from_SP_ID,gsw_cv.Sstar_from_SP);

netcdf.putVar(ncid,CT_chck_cast_ID,gsw_cv.CT_chck_cast);
netcdf.putVar(ncid,CT_from_t_ID,gsw_cv.CT_from_t);

%% other conversions between temperatures, salinities, pressure and height

netcdf.putVar(ncid,deltaSA_from_SP_ID,gsw_cv.deltaSA_from_SP);
netcdf.putVar(ncid,SA_SA_Sstar_from_SP_ID,gsw_cv.SA_SA_Sstar_from_SP);
netcdf.putVar(ncid,Sstar_SA_Sstar_from_SP_ID,gsw_cv.Sstar_SA_Sstar_from_SP);
netcdf.putVar(ncid,SR_from_SP_ID,gsw_cv.SR_from_SP);
netcdf.putVar(ncid,SP_from_SR_ID,gsw_cv.SP_from_SR);
netcdf.putVar(ncid,SP_from_SA_ID,gsw_cv.SP_from_SA);
netcdf.putVar(ncid,Sstar_from_SA_ID,gsw_cv.Sstar_from_SA);
netcdf.putVar(ncid,SA_from_Sstar_ID,gsw_cv.SA_from_Sstar);
netcdf.putVar(ncid,SP_from_Sstar_ID,gsw_cv.SP_from_Sstar);

netcdf.putVar(ncid,t_from_CT_ID,gsw_cv.t_from_CT);
netcdf.putVar(ncid,pt_from_t_ID,gsw_cv.pt_from_t);
netcdf.putVar(ncid,pt0_from_t_ID,gsw_cv.pt0_from_t);
netcdf.putVar(ncid,pt_from_CT_ID,gsw_cv.pt_from_CT);
netcdf.putVar(ncid,CT_from_pt_ID,gsw_cv.CT_from_pt);
netcdf.putVar(ncid,pot_enthalpy_from_pt_ID,gsw_cv.pot_enthalpy_from_pt);

netcdf.putVar(ncid,t90_from_t68_ID,gsw_cv.t90_from_t68);
netcdf.putVar(ncid,t90_from_t48_ID,gsw_cv.t90_from_t48);

netcdf.putVar(ncid,z_from_p_ID,gsw_cv.z_from_p);
netcdf.putVar(ncid,p_from_z_ID,gsw_cv.p_from_z);

netcdf.putVar(ncid,depth_from_z_ID,gsw_cv.depth_from_z);
netcdf.putVar(ncid,z_from_depth_ID,gsw_cv.z_from_depth);

netcdf.putVar(ncid,Abs_Pressure_from_p_ID,gsw_cv.Abs_Pressure_from_p);
netcdf.putVar(ncid,p_from_Abs_Pressure_ID,gsw_cv.p_from_Abs_Pressure);

netcdf.putVar(ncid,entropy_from_CT_ID,gsw_cv.entropy_from_CT);
netcdf.putVar(ncid,CT_from_entropy_ID,gsw_cv.CT_from_entropy);

netcdf.putVar(ncid,entropy_from_pt_ID,gsw_cv.entropy_from_pt);
netcdf.putVar(ncid,pt_from_entropy_ID,gsw_cv.pt_from_entropy);

netcdf.putVar(ncid,entropy_from_t_ID,gsw_cv.entropy_from_t);
netcdf.putVar(ncid,t_from_entropy_ID,gsw_cv.t_from_entropy);

netcdf.putVar(ncid,adiabatic_lapse_rate_from_CT_ID,gsw_cv.adiabatic_lapse_rate_from_CT);
netcdf.putVar(ncid,adiabatic_lapse_rate_from_t_ID,gsw_cv.adiabatic_lapse_rate_from_t);

netcdf.putVar(ncid,molality_from_SA_ID,gsw_cv.molality_from_SA);
netcdf.putVar(ncid,ionic_strength_from_SA_ID,gsw_cv.ionic_strength_from_SA);

%% specific volume, density and enthalpy

netcdf.putVar(ncid,specvol_ID,gsw_cv.specvol);
netcdf.putVar(ncid,alpha_ID,gsw_cv.alpha);
netcdf.putVar(ncid,beta_ID,gsw_cv.beta);
netcdf.putVar(ncid,alpha_on_beta_ID,gsw_cv.alpha_on_beta);

netcdf.putVar(ncid,v_vab_ID,gsw_cv.v_vab);
netcdf.putVar(ncid,alpha_vab_ID,gsw_cv.alpha_vab);
netcdf.putVar(ncid,beta_vab_ID,gsw_cv.beta_vab);

netcdf.putVar(ncid,v_SA_ID,gsw_cv.v_SA);
netcdf.putVar(ncid,v_CT_ID,gsw_cv.v_CT);
netcdf.putVar(ncid,v_P_ID,gsw_cv.v_P);

netcdf.putVar(ncid,v_SA_SA_ID,gsw_cv.v_SA_SA);
netcdf.putVar(ncid,v_SA_CT_ID,gsw_cv.v_SA_CT);
netcdf.putVar(ncid,v_CT_CT_ID,gsw_cv.v_CT_CT);
netcdf.putVar(ncid,v_SA_P_ID,gsw_cv.v_SA_P);
netcdf.putVar(ncid,v_CT_P_ID,gsw_cv.v_CT_P);

netcdf.putVar(ncid,v_SA_wrt_h_ID,gsw_cv.v_SA_wrt_h);
netcdf.putVar(ncid,v_h_ID,gsw_cv.v_h);

netcdf.putVar(ncid,v_SA_SA_wrt_h_ID,gsw_cv.v_SA_SA_wrt_h);
netcdf.putVar(ncid,v_SA_h_ID,gsw_cv.v_SA_h);
netcdf.putVar(ncid,v_h_h_ID,gsw_cv.v_h_h);

netcdf.putVar(ncid,specvol_anom_ID,gsw_cv.specvol_anom);
netcdf.putVar(ncid,specvol_anom_standard_ID,gsw_cv.specvol_anom_standard);

netcdf.putVar(ncid,rho_ID,gsw_cv.rho);

netcdf.putVar(ncid,rho_rab_ID,gsw_cv.rho_rab);
netcdf.putVar(ncid,alpha_rab_ID,gsw_cv.alpha_rab);
netcdf.putVar(ncid,beta_rab_ID,gsw_cv.beta_rab);

netcdf.putVar(ncid,rho_SA_ID,gsw_cv.rho_SA);
netcdf.putVar(ncid,rho_CT_ID,gsw_cv.rho_CT);
netcdf.putVar(ncid,rho_P_ID,gsw_cv.rho_P);

netcdf.putVar(ncid,rho_SA_SA_ID,gsw_cv.rho_SA_SA);
netcdf.putVar(ncid,rho_SA_CT_ID,gsw_cv.rho_SA_CT);
netcdf.putVar(ncid,rho_CT_CT_ID,gsw_cv.rho_CT_CT);
netcdf.putVar(ncid,rho_SA_P_ID,gsw_cv.rho_SA_P);
netcdf.putVar(ncid,rho_CT_P_ID,gsw_cv.rho_CT_P);

netcdf.putVar(ncid,rho_SA_wrt_h_ID,gsw_cv.rho_SA_wrt_h);
netcdf.putVar(ncid,rho_h_ID,gsw_cv.rho_h);

netcdf.putVar(ncid,rho_SA_SA_wrt_h_ID,gsw_cv.rho_SA_SA_wrt_h);
netcdf.putVar(ncid,rho_SA_h_ID,gsw_cv.rho_SA_h);
netcdf.putVar(ncid,rho_h_h_ID,gsw_cv.rho_h_h);

netcdf.putVar(ncid,sigma0_ID,gsw_cv.sigma0);
netcdf.putVar(ncid,sigma1_ID,gsw_cv.sigma1);
netcdf.putVar(ncid,sigma2_ID,gsw_cv.sigma2);
netcdf.putVar(ncid,sigma3_ID,gsw_cv.sigma3);
netcdf.putVar(ncid,sigma4_ID,gsw_cv.sigma4);

netcdf.putVar(ncid,cabbeling_ID,gsw_cv.cabbeling);
netcdf.putVar(ncid,thermobaric_ID,gsw_cv.thermobaric);

netcdf.putVar(ncid,enthalpy_ID,gsw_cv.enthalpy);
netcdf.putVar(ncid,enthalpy_diff_ID,gsw_cv.enthalpy_diff);
netcdf.putVar(ncid,dynamic_enthalpy_ID,gsw_cv.dynamic_enthalpy);

netcdf.putVar(ncid,h_SA_ID,gsw_cv.h_SA);
netcdf.putVar(ncid,h_CT_ID,gsw_cv.h_CT);

netcdf.putVar(ncid,h_SA_SA_ID,gsw_cv.h_SA_SA);
netcdf.putVar(ncid,h_SA_CT_ID,gsw_cv.h_SA_CT);
netcdf.putVar(ncid,h_CT_CT_ID,gsw_cv.h_CT_CT);

netcdf.putVar(ncid,sound_speed_ID,gsw_cv.sound_speed);
netcdf.putVar(ncid,kappa_ID,gsw_cv.kappa);
netcdf.putVar(ncid,internal_energy_ID,gsw_cv.internal_energy);

netcdf.putVar(ncid,u_SA_ID,gsw_cv.u_SA);
netcdf.putVar(ncid,u_CT_ID,gsw_cv.u_CT);
netcdf.putVar(ncid,u_P_ID,gsw_cv.u_P);

netcdf.putVar(ncid,u_SA_SA_ID,gsw_cv.u_SA_SA);
netcdf.putVar(ncid,u_SA_CT_ID,gsw_cv.u_SA_CT);
netcdf.putVar(ncid,u_CT_CT_ID,gsw_cv.u_CT_CT);
netcdf.putVar(ncid,u_SA_P_ID,gsw_cv.u_SA_P);
netcdf.putVar(ncid,u_CT_P_ID,gsw_cv.u_CT_P);

netcdf.putVar(ncid,CT_from_enthalpy_ID,gsw_cv.CT_from_enthalpy);
netcdf.putVar(ncid,SA_from_rho_ID,gsw_cv.SA_from_rho);
netcdf.putVar(ncid,CT_from_rho_ID,gsw_cv.CT_from_rho);
netcdf.putVar(ncid,CT_maxdensity_ID,gsw_cv.CT_maxdensity);

%% vertical stability

netcdf.putVar(ncid,Tu_ID,gsw_cv.Tu);
netcdf.putVar(ncid,Rsubrho_ID,gsw_cv.Rsubrho);
netcdf.putVar(ncid,p_mid_TuRsr_ID,gsw_cv.p_mid_TuRsr);

netcdf.putVar(ncid,n2_ID,gsw_cv.n2);
netcdf.putVar(ncid,p_mid_n2_ID,gsw_cv.p_mid_n2);

netcdf.putVar(ncid,n2min_ID,gsw_cv.n2min);
netcdf.putVar(ncid,n2min_pmid_ID,gsw_cv.n2min_pmid);
netcdf.putVar(ncid,n2min_specvol_ID,gsw_cv.n2min_specvol);
netcdf.putVar(ncid,n2min_alpha_ID,gsw_cv.n2min_alpha);
netcdf.putVar(ncid,n2min_beta_ID,gsw_cv.n2min_beta);
netcdf.putVar(ncid,n2min_dsa_ID,gsw_cv.n2min_dsa);
netcdf.putVar(ncid,n2min_dct_ID,gsw_cv.n2min_dct);
netcdf.putVar(ncid,n2min_dp_ID,gsw_cv.n2min_dp);

netcdf.putVar(ncid,mlp_ID,gsw_cv.mlp);

netcdf.putVar(ncid,n2_lowerlimit_ID,gsw_cv.n2_lowerlimit);

netcdf.putVar(ncid,IPVfN2_ID,gsw_cv.IPVfN2);
netcdf.putVar(ncid,p_mid_IPVfN2_ID,gsw_cv.p_mid_IPVfN2);

%% geostrophic streamfunctions and acoustic travel time

netcdf.putVar(ncid,geo_strf_dyn_height_ID,gsw_cv.geo_strf_dyn_height);

netcdf.putVar(ncid,geo_strf_dyn_height_pc_ID,gsw_cv.geo_strf_dyn_height_pc);
netcdf.putVar(ncid,geo_strf_dyn_height_pc_p_mid_ID,gsw_cv.geo_strf_dyn_height_pc_p_mid);

netcdf.putVar(ncid,geo_strf_isopycnal_ID,gsw_cv.geo_strf_isopycnal);

netcdf.putVar(ncid,geo_strf_isopycnal_pc_ID,gsw_cv.geo_strf_isopycnal_pc);
netcdf.putVar(ncid,geo_strf_isopycnal_pc_p_mid_ID,gsw_cv.geo_strf_isopycnal_pc_p_mid);

netcdf.putVar(ncid,geo_strf_Montgomery_ID,gsw_cv.geo_strf_Montgomery);
netcdf.putVar(ncid,geo_strf_Cunningham_ID,gsw_cv.geo_strf_Cunningham);

netcdf.putVar(ncid,geo_strf_steric_height_ID,gsw_cv.geo_strf_steric_height);

netcdf.putVar(ncid,geo_strf_PISH_ID,gsw_cv.geo_strf_PISH);

netcdf.putVar(ncid,travel_time_ID,gsw_cv.travel_time);

%% Geostrophic velocity

netcdf.putVar(ncid,geo_strf_velocity_ID,gsw_cv.geo_strf_velocity);
netcdf.putVar(ncid,geo_strf_velocity_mid_lat_ID,gsw_cv.geo_strf_velocity_mid_lat);
netcdf.putVar(ncid,geo_strf_velocity_mid_long_ID,gsw_cv.geo_strf_velocity_mid_long);

%% neutral versus isopycnal slopes and ratios
netcdf.putVar(ncid,isopycnal_slope_ratio_ID,gsw_cv.isopycnal_slope_ratio);

netcdf.putVar(ncid,G_CT_ID,gsw_cv.G_CT);
netcdf.putVar(ncid,p_mid_G_CT_ID,gsw_cv.p_mid_G_CT);

netcdf.putVar(ncid,ntpptCT_ID,gsw_cv.ntpptCT);

%% derivatives of enthalpy, entropy, CT and pt

netcdf.putVar(ncid,CT_SA_ID,gsw_cv.CT_SA);
netcdf.putVar(ncid,CT_pt_ID,gsw_cv.CT_pt);

netcdf.putVar(ncid,CT_SA_SA_ID,gsw_cv.CT_SA_SA);
netcdf.putVar(ncid,CT_SA_pt_ID,gsw_cv.CT_SA_pt);
netcdf.putVar(ncid,CT_pt_pt_ID,gsw_cv.CT_pt_pt);

netcdf.putVar(ncid,eta_SA_ID,gsw_cv.eta_SA);
netcdf.putVar(ncid,eta_CT_ID,gsw_cv.eta_CT);

netcdf.putVar(ncid,eta_SA_SA_ID,gsw_cv.eta_SA_SA);
netcdf.putVar(ncid,eta_SA_CT_ID,gsw_cv.eta_SA_CT);
netcdf.putVar(ncid,eta_CT_CT_ID,gsw_cv.eta_CT_CT);

netcdf.putVar(ncid,pt_SA_ID,gsw_cv.pt_SA);
netcdf.putVar(ncid,pt_CT_ID,gsw_cv.pt_CT);

netcdf.putVar(ncid,pt_SA_SA_ID,gsw_cv.pt_SA_SA);
netcdf.putVar(ncid,pt_SA_CT_ID,gsw_cv.pt_SA_CT);
netcdf.putVar(ncid,pt_CT_CT_ID,gsw_cv.pt_CT_CT);

%% seawater properties at freezing temperatures

netcdf.putVar(ncid,CT_freezing_ID,gsw_cv.CT_freezing);
netcdf.putVar(ncid,CT_freezing_poly_ID,gsw_cv.CT_freezing_poly);
netcdf.putVar(ncid,t_freezing_ID,gsw_cv.t_freezing);
netcdf.putVar(ncid,t_freezing_poly_ID,gsw_cv.t_freezing_poly);
netcdf.putVar(ncid,pot_enthalpy_ice_freezing_ID,gsw_cv.pot_enthalpy_ice_freezing);
netcdf.putVar(ncid,pot_enthalpy_ice_freezing_poly_ID,gsw_cv.pot_enthalpy_ice_freezing_poly);

netcdf.putVar(ncid,SA_freezing_from_CT_ID,gsw_cv.SA_freezing_from_CT);
netcdf.putVar(ncid,SA_freezing_from_CT_poly_ID,gsw_cv.SA_freezing_from_CT_poly);
netcdf.putVar(ncid,SA_freezing_from_t_ID,gsw_cv.SA_freezing_from_t);
netcdf.putVar(ncid,SA_freezing_from_t_poly_ID,gsw_cv.SA_freezing_from_t_poly);

netcdf.putVar(ncid,pressure_freezing_CT_ID,gsw_cv.pressure_freezing_CT);

netcdf.putVar(ncid,CTfreezing_SA_ID,gsw_cv.CTfreezing_SA);
netcdf.putVar(ncid,CTfreezing_P_ID,gsw_cv.CTfreezing_P);

netcdf.putVar(ncid,CTfreezing_SA_poly_ID,gsw_cv.CTfreezing_SA_poly);
netcdf.putVar(ncid,CTfreezing_P_poly_ID,gsw_cv.CTfreezing_P_poly);

netcdf.putVar(ncid,tfreezing_SA_ID,gsw_cv.tfreezing_SA);
netcdf.putVar(ncid,tfreezing_P_ID,gsw_cv.tfreezing_P);

netcdf.putVar(ncid,tfreezing_SA_poly_ID,gsw_cv.tfreezing_SA_poly);
netcdf.putVar(ncid,tfreezing_P_poly_ID,gsw_cv.tfreezing_P_poly);

netcdf.putVar(ncid,pot_enthalpy_ice_freezing_SA_ID,gsw_cv.pot_enthalpy_ice_freezing_SA);
netcdf.putVar(ncid,pot_enthalpy_ice_freezing_P_ID,gsw_cv.pot_enthalpy_ice_freezing_P);

netcdf.putVar(ncid,pot_enthalpy_ice_freezing_SA_poly_ID,gsw_cv.pot_enthalpy_ice_freezing_SA_poly);
netcdf.putVar(ncid,pot_enthalpy_ice_freezing_P_poly_ID,gsw_cv.pot_enthalpy_ice_freezing_P_poly);

netcdf.putVar(ncid,latentheat_melting_ID,gsw_cv.latentheat_melting);

%%  thermodynamic interaction between ice and seawater

netcdf.putVar(ncid,melting_ice_SA_CT_ratio_ID,gsw_cv.melting_ice_SA_CT_ratio);
netcdf.putVar(ncid,melting_ice_SA_CT_ratio_poly_ID,gsw_cv.melting_ice_SA_CT_ratio_poly);

netcdf.putVar(ncid,melting_ice_equilibrium_SA_CT_ratio_ID,gsw_cv.melting_ice_equilibrium_SA_CT_ratio);
netcdf.putVar(ncid,melting_ice_equilibrium_SA_CT_ratio_poly_ID,gsw_cv.melting_ice_equilibrium_SA_CT_ratio_poly);

netcdf.putVar(ncid,melting_ice_into_seawater_SA_final_ID,gsw_cv.melting_ice_into_seawater_SA_final);
netcdf.putVar(ncid,melting_ice_into_seawater_CT_final_ID,gsw_cv.melting_ice_into_seawater_CT_final);

netcdf.putVar(ncid,ice_fraction_to_freeze_seawater_SA_freeze_ID,gsw_cv.ice_fraction_to_freeze_seawater_SA_freeze);
netcdf.putVar(ncid,ice_fraction_to_freeze_seawater_CT_freeze_ID,gsw_cv.ice_fraction_to_freeze_seawater_CT_freeze);
netcdf.putVar(ncid,ice_fraction_to_freeze_seawater_w_Ih_ID,gsw_cv.ice_fraction_to_freeze_seawater_w_Ih);

netcdf.putVar(ncid,dSA_dCT_frazil_ID,gsw_cv.dSA_dCT_frazil);
netcdf.putVar(ncid,dSA_dP_frazil_ID,gsw_cv.dSA_dP_frazil);
netcdf.putVar(ncid,dCT_dP_frazil_ID,gsw_cv.dCT_dP_frazil);

netcdf.putVar(ncid,dSA_dCT_frazil_poly_ID,gsw_cv.dSA_dCT_frazil_poly);
netcdf.putVar(ncid,dSA_dP_frazil_poly_ID,gsw_cv.dSA_dP_frazil_poly);
netcdf.putVar(ncid,dCT_dP_frazil_poly_ID,gsw_cv.dCT_dP_frazil_poly);

netcdf.putVar(ncid,frazil_properties_potential_SA_final_ID,gsw_cv.frazil_properties_potential_SA_final);
netcdf.putVar(ncid,frazil_properties_potential_CT_final_ID,gsw_cv.frazil_properties_potential_CT_final);
netcdf.putVar(ncid,frazil_properties_potential_w_Ih_final_ID,gsw_cv.frazil_properties_potential_w_Ih_final);

netcdf.putVar(ncid,frazil_properties_potential_poly_SA_final_ID,gsw_cv.frazil_properties_potential_poly_SA_final);
netcdf.putVar(ncid,frazil_properties_potential_poly_CT_final_ID,gsw_cv.frazil_properties_potential_poly_CT_final);
netcdf.putVar(ncid,frazil_properties_potential_poly_w_Ih_final_ID,gsw_cv.frazil_properties_potential_poly_w_Ih_final);

netcdf.putVar(ncid,frazil_properties_SA_final_ID,gsw_cv.frazil_properties_SA_final);
netcdf.putVar(ncid,frazil_properties_CT_final_ID,gsw_cv.frazil_properties_CT_final);
netcdf.putVar(ncid,frazil_properties_w_Ih_final_ID,gsw_cv.frazil_properties_w_Ih_final);

%%  thermodynamic interaction between sea ice and seawater

netcdf.putVar(ncid,melting_seaice_SA_CT_ratio_ID,gsw_cv.melting_seaice_SA_CT_ratio);
netcdf.putVar(ncid,melting_seaice_SA_CT_ratio_poly_ID,gsw_cv.melting_seaice_SA_CT_ratio_poly);

netcdf.putVar(ncid,melting_seaice_equilibrium_SA_CT_ratio_ID,gsw_cv.melting_seaice_equilibrium_SA_CT_ratio);
netcdf.putVar(ncid,melting_seaice_equilibrium_SA_CT_ratio_poly_ID,gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_poly);

netcdf.putVar(ncid,melting_seaice_into_seawater_SA_final_ID,gsw_cv.melting_seaice_into_seawater_SA_final);
netcdf.putVar(ncid,melting_seaice_into_seawater_CT_final_ID,gsw_cv.melting_seaice_into_seawater_CT_final);

netcdf.putVar(ncid,seaice_fraction_to_freeze_seawater_SA_freeze_ID,gsw_cv.seaice_fraction_to_freeze_seawater_SA_freeze);
netcdf.putVar(ncid,seaice_fraction_to_freeze_seawater_CT_freeze_ID,gsw_cv.seaice_fraction_to_freeze_seawater_CT_freeze);
netcdf.putVar(ncid,seaice_fraction_to_freeze_seawater_w_Ih_ID,gsw_cv.seaice_fraction_to_freeze_seawater_w_Ih);

%% themodynamic properties of ice Ih

netcdf.putVar(ncid,specvol_ice_ID,gsw_cv.specvol_ice);

netcdf.putVar(ncid,alpha_wrt_t_ice_ID,gsw_cv.alpha_wrt_t_ice);

netcdf.putVar(ncid,rho_ice_ID,gsw_cv.rho_ice);

netcdf.putVar(ncid,pressure_coefficient_ice_ID,gsw_cv.pressure_coefficient_ice);
netcdf.putVar(ncid,sound_speed_ice_ID,gsw_cv.sound_speed_ice);
netcdf.putVar(ncid,kappa_ice_ID,gsw_cv.kappa_ice);
netcdf.putVar(ncid,kappa_const_t_ice_ID,gsw_cv.kappa_const_t_ice);
netcdf.putVar(ncid,internal_energy_ice_ID,gsw_cv.internal_energy_ice);

netcdf.putVar(ncid,enthalpy_ice_ID,gsw_cv.enthalpy_ice);

netcdf.putVar(ncid,entropy_ice_ID,gsw_cv.entropy_ice);

netcdf.putVar(ncid,cp_ice_ID,gsw_cv.cp_ice);
netcdf.putVar(ncid,chem_potential_water_ice_ID,gsw_cv.chem_potential_water_ice);
netcdf.putVar(ncid,Helmholtz_energy_ice_ID,gsw_cv.Helmholtz_energy_ice);

netcdf.putVar(ncid,adiabatic_lapse_rate_ice_ID,gsw_cv.adiabatic_lapse_rate_ice);

netcdf.putVar(ncid,pt0_from_t_ice_ID,gsw_cv.pt0_from_t_ice);
netcdf.putVar(ncid,pt_from_t_ice_ID,gsw_cv.pt_from_t_ice);
netcdf.putVar(ncid,t_from_pt0_ice_ID,gsw_cv.t_from_pt0_ice);
netcdf.putVar(ncid,t_from_rho_ice_ID,gsw_cv.t_from_rho_ice);

netcdf.putVar(ncid,pot_enthalpy_from_pt_ice_ID,gsw_cv.pot_enthalpy_from_pt_ice);
netcdf.putVar(ncid,pt_from_pot_enthalpy_ice_ID,gsw_cv.pt_from_pot_enthalpy_ice);

netcdf.putVar(ncid,pot_enthalpy_from_pt_ice_poly_ID,gsw_cv.pot_enthalpy_from_pt_ice_poly);
netcdf.putVar(ncid,pt_from_pot_enthalpy_ice_poly_ID,gsw_cv.pt_from_pot_enthalpy_ice_poly);

netcdf.putVar(ncid,pot_enthalpy_from_specvol_ice_ID,gsw_cv.pot_enthalpy_from_specvol_ice);
netcdf.putVar(ncid,specvol_from_pot_enthalpy_ice_ID,gsw_cv.specvol_from_pot_enthalpy_ice);

netcdf.putVar(ncid,pot_enthalpy_from_specvol_ice_poly_ID,gsw_cv.pot_enthalpy_from_specvol_ice_poly);
netcdf.putVar(ncid,specvol_from_pot_enthalpy_ice_poly_ID,gsw_cv.specvol_from_pot_enthalpy_ice_poly);

%% isobaric evaporation enthalpy 

netcdf.putVar(ncid,latentheat_evap_CT_ID,gsw_cv.latentheat_evap_CT);
netcdf.putVar(ncid,latentheat_evap_t_ID,gsw_cv.latentheat_evap_t);

%% spiciness

netcdf.putVar(ncid,spiciness0_ID,gsw_cv.spiciness0);
netcdf.putVar(ncid,spiciness1_ID,gsw_cv.spiciness1);
netcdf.putVar(ncid,spiciness2_ID,gsw_cv.spiciness2);

%% planet earth properties

netcdf.putVar(ncid,f_ID,gsw_cv.f);
netcdf.putVar(ncid,grav_ID,gsw_cv.grav);
netcdf.putVar(ncid,distance_ID,gsw_cv.distance);

%% TEOS-10 constants

netcdf.putVar(ncid,T0_ID,gsw_cv.T0);
netcdf.putVar(ncid,P0_ID,gsw_cv.P0);
netcdf.putVar(ncid,SSO_ID,gsw_cv.SSO);
netcdf.putVar(ncid,uPS_ID,gsw_cv.uPS);
netcdf.putVar(ncid,cp0_ID,gsw_cv.cp0);
netcdf.putVar(ncid,C3515_ID,gsw_cv.C3515);
netcdf.putVar(ncid,SonCl_ID,gsw_cv.SonCl);
netcdf.putVar(ncid,valence_factor_ID,gsw_cv.valence_factor);
netcdf.putVar(ncid,atomic_weight_ID,gsw_cv.atomic_weight);

%% dissolved gasses
netcdf.putVar(ncid,Arsol_ID,gsw_cv.Arsol);
netcdf.putVar(ncid,Arsol_SP_pt_ID,gsw_cv.Arsol_SP_pt);

netcdf.putVar(ncid,Hesol_ID,gsw_cv.Hesol);
netcdf.putVar(ncid,Hesol_SP_pt_ID,gsw_cv.Hesol_SP_pt);

netcdf.putVar(ncid,Krsol_ID,gsw_cv.Krsol);
netcdf.putVar(ncid,Krsol_SP_pt_ID,gsw_cv.Krsol_SP_pt);

netcdf.putVar(ncid,N2Osol_ID,gsw_cv.N2Osol);
netcdf.putVar(ncid,N2Osol_SP_pt_ID,gsw_cv.N2Osol_SP_pt);

netcdf.putVar(ncid,N2sol_ID,gsw_cv.N2sol);
netcdf.putVar(ncid,N2sol_SP_pt_ID,gsw_cv.N2sol_SP_pt);

netcdf.putVar(ncid,Nesol_ID,gsw_cv.Nesol);
netcdf.putVar(ncid,Nesol_SP_pt_ID,gsw_cv.Nesol_SP_pt);

netcdf.putVar(ncid,O2sol_ID,gsw_cv.O2sol);
netcdf.putVar(ncid,O2sol_SP_pt_ID,gsw_cv.O2sol_SP_pt);

%% density and enthalpy in terms of CT, derived from the exact Gibbs function

netcdf.putVar(ncid,specvol_CT_exact_ID,gsw_cv.specvol_CT_exact);

netcdf.putVar(ncid,alpha_CT_exact_ID,gsw_cv.alpha_CT_exact);
netcdf.putVar(ncid,beta_CT_exact_ID,gsw_cv.beta_CT_exact);
netcdf.putVar(ncid,alpha_on_beta_CT_exact_ID,gsw_cv.alpha_on_beta_CT_exact);

netcdf.putVar(ncid,v_vab_CT_exact_ID,gsw_cv.v_vab_CT_exact);
netcdf.putVar(ncid,alpha_vab_CT_exact_ID,gsw_cv.alpha_vab_CT_exact);
netcdf.putVar(ncid,beta_vab_CT_exact_ID,gsw_cv.beta_vab_CT_exact);

netcdf.putVar(ncid,v_SA_CT_exact_ID,gsw_cv.v_SA_CT_exact);
netcdf.putVar(ncid,v_CT_CT_exact_ID,gsw_cv.v_CT_CT_exact);
netcdf.putVar(ncid,v_P_CT_exact_ID,gsw_cv.v_P_CT_exact);

netcdf.putVar(ncid,v_SA_SA_CT_exact_ID,gsw_cv.v_SA_SA_CT_exact);
netcdf.putVar(ncid,v_SA_CT_CT_exact_ID,gsw_cv.v_SA_CT_CT_exact);
netcdf.putVar(ncid,v_CT_CT_CT_exact_ID,gsw_cv.v_CT_CT_CT_exact);
netcdf.putVar(ncid,v_SA_P_CT_exact_ID,gsw_cv.v_SA_P_CT_exact);
netcdf.putVar(ncid,v_CT_P_CT_exact_ID,gsw_cv.v_CT_P_CT_exact);

netcdf.putVar(ncid,v_SA_wrt_h_CT_exact_ID,gsw_cv.v_SA_wrt_h_CT_exact);
netcdf.putVar(ncid,v_h_CT_exact_ID,gsw_cv.v_h_CT_exact);

netcdf.putVar(ncid,v_SA_SA_wrt_h_CT_exact_ID,gsw_cv.v_SA_SA_wrt_h_CT_exact);
netcdf.putVar(ncid,v_SA_h_CT_exact_ID,gsw_cv.v_SA_h_CT_exact);
netcdf.putVar(ncid,v_h_h_CT_exact_ID,gsw_cv.v_h_h_CT_exact);

netcdf.putVar(ncid,specvol_anom_CT_exact_ID,gsw_cv.specvol_anom_CT_exact);
netcdf.putVar(ncid,specvol_anom_standard_CT_exact_ID,gsw_cv.specvol_anom_standard_CT_exact);

netcdf.putVar(ncid,rho_CT_exact_ID,gsw_cv.rho_CT_exact);

netcdf.putVar(ncid,rho_rab_CT_exact_ID,gsw_cv.rho_rab_CT_exact);
netcdf.putVar(ncid,alpha_rab_CT_exact_ID,gsw_cv.alpha_rab_CT_exact);
netcdf.putVar(ncid,beta_rab_CT_exact_ID,gsw_cv.beta_rab_CT_exact);

netcdf.putVar(ncid,rho_SA_CT_exact_ID,gsw_cv.rho_SA_CT_exact);
netcdf.putVar(ncid,rho_CT_CT_exact_ID,gsw_cv.rho_CT_CT_exact);
netcdf.putVar(ncid,rho_P_CT_exact_ID,gsw_cv.rho_P_CT_exact);

netcdf.putVar(ncid,rho_SA_SA_CT_exact_ID,gsw_cv.rho_SA_SA_CT_exact);
netcdf.putVar(ncid,rho_SA_CT_CT_exact_ID,gsw_cv.rho_SA_CT_CT_exact);
netcdf.putVar(ncid,rho_CT_CT_CT_exact_ID,gsw_cv.rho_CT_CT_CT_exact);
netcdf.putVar(ncid,rho_SA_P_CT_exact_ID,gsw_cv.rho_SA_P_CT_exact);
netcdf.putVar(ncid,rho_CT_P_CT_exact_ID,gsw_cv.rho_CT_P_CT_exact);

netcdf.putVar(ncid,rho_SA_wrt_h_CT_exact_ID,gsw_cv.rho_SA_wrt_h_CT_exact);
netcdf.putVar(ncid,rho_h_CT_exact_ID,gsw_cv.rho_h_CT_exact);

netcdf.putVar(ncid,rho_SA_SA_wrt_h_CT_exact_ID,gsw_cv.rho_SA_SA_wrt_h_CT_exact);
netcdf.putVar(ncid,rho_SA_h_CT_exact_ID,gsw_cv.rho_SA_h_CT_exact);
netcdf.putVar(ncid,rho_h_h_CT_exact_ID,gsw_cv.rho_h_h_CT_exact);

netcdf.putVar(ncid,sigma0_CT_exact_ID,gsw_cv.sigma0_CT_exact);
netcdf.putVar(ncid,sigma1_CT_exact_ID,gsw_cv.sigma1_CT_exact);
netcdf.putVar(ncid,sigma2_CT_exact_ID,gsw_cv.sigma2_CT_exact);
netcdf.putVar(ncid,sigma3_CT_exact_ID,gsw_cv.sigma3_CT_exact);
netcdf.putVar(ncid,sigma4_CT_exact_ID,gsw_cv.sigma4_CT_exact);

netcdf.putVar(ncid,cabbeling_CT_exact_ID,gsw_cv.cabbeling_CT_exact);
netcdf.putVar(ncid,thermobaric_CT_exact_ID,gsw_cv.thermobaric_CT_exact);

netcdf.putVar(ncid,enthalpy_CT_exact_ID,gsw_cv.enthalpy_CT_exact);
netcdf.putVar(ncid,enthalpy_diff_CT_exact_ID,gsw_cv.enthalpy_diff_CT_exact);
netcdf.putVar(ncid,dynamic_enthalpy_CT_exact_ID,gsw_cv.dynamic_enthalpy_CT_exact);

netcdf.putVar(ncid,h_SA_CT_exact_ID,gsw_cv.h_SA_CT_exact);
netcdf.putVar(ncid,h_CT_CT_exact_ID,gsw_cv.h_CT_CT_exact);

netcdf.putVar(ncid,h_SA_SA_CT_exact_ID,gsw_cv.h_SA_SA_CT_exact);
netcdf.putVar(ncid,h_SA_CT_CT_exact_ID,gsw_cv.h_SA_CT_CT_exact);
netcdf.putVar(ncid,h_CT_CT_CT_exact_ID,gsw_cv.h_CT_CT_CT_exact);

netcdf.putVar(ncid,sound_speed_CT_exact_ID,gsw_cv.sound_speed_CT_exact);
netcdf.putVar(ncid,kappa_CT_exact_ID,gsw_cv.kappa_CT_exact);

netcdf.putVar(ncid,internal_energy_CT_exact_ID,gsw_cv.internal_energy_CT_exact);

netcdf.putVar(ncid,u_SA_CT_exact_ID,gsw_cv.u_SA_CT_exact);
netcdf.putVar(ncid,u_CT_CT_exact_ID,gsw_cv.u_CT_CT_exact);
netcdf.putVar(ncid,u_P_CT_exact_ID,gsw_cv.u_P_CT_exact);

netcdf.putVar(ncid,u_SA_SA_CT_exact_ID,gsw_cv.u_SA_SA_CT_exact);
netcdf.putVar(ncid,u_SA_CT_CT_exact_ID,gsw_cv.u_SA_CT_CT_exact);
netcdf.putVar(ncid,u_CT_CT_CT_exact_ID,gsw_cv.u_CT_CT_CT_exact);
netcdf.putVar(ncid,u_SA_P_CT_exact_ID,gsw_cv.u_SA_P_CT_exact);
netcdf.putVar(ncid,u_CT_P_CT_exact_ID,gsw_cv.u_CT_P_CT_exact);

netcdf.putVar(ncid,CT_from_enthalpy_ID,gsw_cv.CT_from_enthalpy);
netcdf.putVar(ncid,SA_from_rho_CT_exact_ID,gsw_cv.SA_from_rho_CT_exact);
netcdf.putVar(ncid,CT_maxdensity_exact_ID,gsw_cv.CT_maxdensity_exact);
netcdf.putVar(ncid,CT_from_rho_exact_ID,gsw_cv.CT_from_rho_exact);

%% Labrortory functions

netcdf.putVar(ncid,rho_t_exact_ID,gsw_cv.rho_t_exact);
netcdf.putVar(ncid,SA_from_rho_t_exact_ID,gsw_cv.SA_from_rho_t_exact);
netcdf.putVar(ncid,deltaSA_from_rho_t_exact_ID,gsw_cv.deltaSA_from_rho_t_exact);

%% basic thermodynamic properties interms of in-situ t, derived from the exact Gibbs function

netcdf.putVar(ncid,specvol_t_exact_ID,gsw_cv.specvol_t_exact);

netcdf.putVar(ncid,alpha_wrt_CT_t_exact_ID,gsw_cv.alpha_wrt_CT_t_exact);
netcdf.putVar(ncid,alpha_wrt_pt_t_exact_ID,gsw_cv.alpha_wrt_pt_t_exact);
netcdf.putVar(ncid,alpha_wrt_t_exact_ID,gsw_cv.alpha_wrt_t_exact);

netcdf.putVar(ncid,beta_const_CT_t_exact_ID,gsw_cv.beta_const_CT_t_exact);
netcdf.putVar(ncid,beta_const_pt_t_exact_ID,gsw_cv.beta_const_pt_t_exact);
netcdf.putVar(ncid,beta_const_t_exact_ID,gsw_cv.beta_const_t_exact);

netcdf.putVar(ncid,specvol_anom_standard_t_exact_ID,gsw_cv.specvol_anom_standard_t_exact);

netcdf.putVar(ncid,pot_rho_t_exact_ID,gsw_cv.pot_rho_t_exact);
netcdf.putVar(ncid,sigma0_pt0_exact_ID,gsw_cv.sigma0_pt0_exact);

netcdf.putVar(ncid,enthalpy_t_exact_ID,gsw_cv.enthalpy_t_exact);
netcdf.putVar(ncid,dynamic_enthalpy_t_exact_ID,gsw_cv.dynamic_enthalpy_t_exact);

netcdf.putVar(ncid,CT_SA_wrt_t_ID,gsw_cv.CT_SA_wrt_t);
netcdf.putVar(ncid,CT_T_wrt_t_ID,gsw_cv.CT_T_wrt_t);
netcdf.putVar(ncid,CT_P_wrt_t_ID,gsw_cv.CT_P_wrt_t);

netcdf.putVar(ncid,h_SA_wrt_t_ID,gsw_cv.h_SA_wrt_t);
netcdf.putVar(ncid,h_T_wrt_t_ID,gsw_cv.h_T_wrt_t);
netcdf.putVar(ncid,h_P_wrt_t_ID,gsw_cv.h_P_wrt_t);

netcdf.putVar(ncid,sound_speed_t_exact_ID,gsw_cv.sound_speed_t_exact);
netcdf.putVar(ncid,kappa_t_exact_ID,gsw_cv.kappa_t_exact);
netcdf.putVar(ncid,kappa_const_t_exact_ID,gsw_cv.kappa_const_t_exact);
netcdf.putVar(ncid,internal_energy_t_exact_ID,gsw_cv.internal_energy_t_exact);

netcdf.putVar(ncid,t_from_rho_exact_ID,gsw_cv.t_from_rho_exact);
netcdf.putVar(ncid,t_maxdensity_exact_ID,gsw_cv.t_maxdensity_exact);

netcdf.putVar(ncid,cp_t_exact_ID,gsw_cv.cp_t_exact);

netcdf.putVar(ncid,isochoric_heat_cap_t_exact_ID,gsw_cv.isochoric_heat_cap_t_exact);

netcdf.putVar(ncid,chem_potential_relative_t_exact_ID,gsw_cv.chem_potential_relative_t_exact);
netcdf.putVar(ncid,chem_potential_water_t_exact_ID,gsw_cv.chem_potential_water_t_exact);
netcdf.putVar(ncid,chem_potential_salt_t_exact_ID,gsw_cv.chem_potential_salt_t_exact);
netcdf.putVar(ncid,t_deriv_chem_potential_water_t_exact_ID,gsw_cv.t_deriv_chem_potential_water_t_exact);

netcdf.putVar(ncid,dilution_coefficient_t_exact_ID,gsw_cv.dilution_coefficient_t_exact);

netcdf.putVar(ncid,Gibbs_energy_t_exact_ID,gsw_cv.Gibbs_energy_t_exact);
netcdf.putVar(ncid,Helmholtz_energy_t_exact_ID,gsw_cv.Helmholtz_energy_t_exact);
netcdf.putVar(ncid,osmotic_coefficient_t_exact_ID,gsw_cv.osmotic_coefficient_t_exact);
netcdf.putVar(ncid,osmotic_pressure_t_exact_ID,gsw_cv.osmotic_pressure_t_exact);

%% Library

netcdf.putVar(ncid,Fdelta_ID,gsw_cv.Fdelta);
netcdf.putVar(ncid,deltaSA_atlas_ID,gsw_cv.deltaSA_atlas);

%We're done, close the netcdf
netcdf.close(ncid)

