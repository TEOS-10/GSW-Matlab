% GSW Oceanographic Toolbox 
% Version 3.06.12 (R2018a) 18-June-2020
%
% documentation set
%  gsw_front_page                                       - front page to the GSW Oceanographic Toolbox
%  gsw_check_functions                                  - checks that all the GSW functions work correctly
%  gsw_demo                                             - demonstrates many GSW functions and features
%  gsw_ver                                              - displays the GSW version number
%  gsw_licence                                          - creative commons licence for the GSW Oceanographic Toolbox
%
% Practical Salinity (SP), PSS-78  
%  gsw_SP_from_C                                        - Practical Salinity from conductivity, C (inc. for SP < 2)
%  gsw_C_from_SP                                        - conductivity, C, from Practical Salinity (inc. for SP < 2)
%  gsw_SP_from_R                                        - Practical Salinity from conductivity ratio, R (inc. for SP < 2)
%  gsw_R_from_SP                                        - conductivity ratio, R, from Practical Salinity (inc. for SP < 2)
%  gsw_SP_salinometer                                   - Practical Salinity from a laboratory salinometer (inc. for SP < 2)
%  gsw_SP_from_SK                                       - Practical Salinity from Knudsen Salinity
%
% Absolute Salinity (SA), Preformed Salinity (Sstar) and Conservative Temperature (CT) 
%  gsw_SA_from_SP                                       - Absolute Salinity from Practical Salinity
%  gsw_Sstar_from_SP                                    - Preformed Salinity from Practical Salinity
%  gsw_CT_from_t                                        - Conservative Temperature from in-situ temperature
%
% Absolute Salinity - Conservative Temperature plotting function 
%  gsw_SA_CT_plot                                       - function to plot Absolute Salinity - Conservative Temperature
%                                                         profiles on the SA-CT diagram, including the freezing line
%                                                         and selected potential density contours
%
% other conversions between temperatures, salinities, pressure and height
%  gsw_deltaSA_from_SP                                  - Absolute Salinity Anomaly from Practical Salinity
%  gsw_SA_Sstar_from_SP                                 - Absolute Salinity & Preformed Salinity from Practical Salinity
%  gsw_SR_from_SP                                       - Reference Salinity from Practical Salinity
%  gsw_SP_from_SR                                       - Practical Salinity from Reference Salinity
%  gsw_SP_from_SA                                       - Practical Salinity from Absolute Salinity
%  gsw_Sstar_from_SA                                    - Preformed Salinity from Absolute Salinity
%  gsw_SA_from_Sstar                                    - Absolute Salinity from Preformed Salinity
%  gsw_SP_from_Sstar                                    - Practical Salinity from Preformed Salinity
%  gsw_pt_from_CT                                       - potential temperature from Conservative Temperature
%  gsw_t_from_CT                                        - in-situ temperature from Conservative Temperature
%  gsw_CT_from_pt                                       - Conservative Temperature from potential temperature
%  gsw_pot_enthalpy_from_pt                             - potential enthalpy from potential temperature
%  gsw_pt_from_t                                        - potential temperature
%  gsw_pt0_from_t                                       - potential temperature with a reference pressure of zero dbar
%  gsw_t_from_pt0                                       - in-situ temperature from potential temperature with p_ref = 0 dbar
%  gsw_t90_from_t48                                     - ITS-90 temperature from IPTS-48 temperature
%  gsw_t90_from_t68                                     - ITS-90 temperature from IPTS-68 temperature
%  gsw_z_from_p                                         - height from pressure
%  gsw_p_from_z                                         - pressure from height
%  gsw_z_from_depth                                     - height from depth
%  gsw_depth_from_z                                     - depth from height
%  gsw_Abs_Pressure_from_p                              - Absolute Pressure,P, from pressure, p
%  gsw_p_from_Abs_Pressure                              - pressure, p, from Absolute Pressure, P
%  gsw_entropy_from_CT                                  - entropy from Conservative Temperature
%  gsw_CT_from_entropy                                  - Conservative Temperature from entropy 
%  gsw_entropy_from_pt                                  - entropy from potential temperature
%  gsw_pt_from_entropy                                  - potential temperature from entropy
%  gsw_entropy_from_t                                   - entropy from in-situ temperature
%  gsw_t_from_entropy                                   - in-situ temperature from entropy
%  gsw_adiabatic_lapse_rate_from_CT                     - adiabatic lapse rate from Conservative Temperature
%  gsw_adiabatic_lapse_rate_from_t                      - adiabatic lapse rate from in-situ temperature
%  gsw_molality_from_SA                                 - molality of seawater
%  gsw_ionic_strength_from_SA                           - ionic strength of seawater
%
% specific volume, density and enthalpy 
%  gsw_specvol                                          - specific volume
%  gsw_alpha                                            - thermal expansion coefficient with respect to CT
%  gsw_beta                                             - saline contraction coefficient at constant CT
%  gsw_alpha_on_beta                                    - alpha divided by beta
%  gsw_specvol_alpha_beta                               - specific volume, thermal expansion & saline contraction coefficients
%  gsw_specvol_first_derivatives                        - first derivatives of specific volume
%  gsw_specvol_second_derivatives                       - second derivatives of specific volume
%  gsw_specvol_first_derivatives_wrt_enthalpy           - first derivatives of specific volume
%  gsw_specvol_second_derivatives_wrt_enthalpy          - second derivatives of specific volume
%  gsw_specvol_anom                                     - specific volume anomaly
%  gsw_specvol_anom_standard                            - specific volume anomaly relative to SSO & 0C
%  gsw_specvol_diff                                     - specific volume difference between two pressures
%  gsw_rho                                              - in-situ density and potential density
%  gsw_rho_alpha_beta                                   - in-situ density, thermal expansion & saline contraction coefficients
%  gsw_rho_first_derivatives                            - first derivatives of density
%  gsw_rho_first_derivatives                            - first derivatives of density
%  gsw_rho_second_derivatives                           - second derivatives of density
%  gsw_rho_first_derivatives_wrt_enthalpy               - first derivatives of density
%  gsw_rho_second_derivatives_wrt_enthalpy              - second derivatives of density
%  gsw_sigma0                                           - sigma0 with reference pressure of 0 dbar
%  gsw_sigma1                                           - sigma1 with reference pressure of 1000 dbar
%  gsw_sigma2                                           - sigma2 with reference pressure of 2000 dbar
%  gsw_sigma3                                           - sigma3 with reference pressure of 3000 dbar
%  gsw_sigma4                                           - sigma4 with reference pressure of 4000 dbar
%  gsw_cabbeling                                        - cabbeling coefficient
%  gsw_thermobaric                                      - thermobaric coefficient
%  gsw_enthalpy                                         - enthalpy
%  gsw_enthalpy_diff                                    - difference of enthalpy between two pressures
%  gsw_dynamic_enthalpy                                 - dynamic enthalpy
%  gsw_enthalpy_first_derivatives                       - first derivatives of enthalpy
%  gsw_enthalpy_second_derivatives                      - second derivatives of enthalpy
%  gsw_sound_speed                                      - sound speed (approximate, with r.m.s. error of 0.067 m/s)
%  gsw_kappa                                            - isentropic compressibility
%  gsw_internal_energy                                  - internal energy
%  gsw_internal_energy_first_derivatives                - first derivatives of internal energy
%  gsw_internal_energy_second_derivatives               - second derivatives of internal energy
%  gsw_CT_from_enthalpy                                 - Conservative Temperature from enthalpy
%  gsw_SA_from_rho                                      - Absolute Salinity from density
%  gsw_CT_from_rho                                      - Conservative Temperature from density
%  gsw_CT_maxdensity                                    - Conservative Temperature of maximum density of seawater
%
% vertical stability and interpolation
%  gsw_Turner_Rsubrho                                   - Turner angle & Rsubrho
%  gsw_Nsquared                                         - buoyancy (Brunt-Vaisala) frequency squared (N^2)
%  gsw_Nsquared_min                                     - minimum buoyancy (Brunt-Vaisala) frequency squared (N^2)
%  gsw_stabilise_SA_const_t	                            - minimally adjust SA to produce a stable water column, keeping in-situ temperature constant
%  gsw_stabilise_SA_CT	                                - minimally adjusts SA & CT to produce a stable water column
%  gsw_mlp	                                            - mixed-layer pressure
%  gsw_Nsquared_lowerlimit                              - specified profile of minimum buoyancy frequency squared
%  gsw_SA_CT_interp                                     - interpolates SA & CT to defined pressures
%  gsw_t_interp                                         - interpolates in-situ temperature to defined pressures
%  gsw_tracer_CT_interp                                 - interpolates a tracer & CT to defined pressures
%  gsw_tracer_interp                                    - interpolates a tracer to defined pressures
%  gsw_IPV_vs_fNsquared_ratio                           - ratio of the vertical gradient of potential density
%                                                         (with reference pressure, p_ref), to the vertical 
%                                                         gradient of locally-referenced potential density
%
% geostrophic streamfunctions, acoustic travel time and geostrophic velocity
%  gsw_geo_strf_dyn_height                              - dynamic height anomaly
%  gsw_geo_strf_dyn_height_pc                           - dynamic height anomaly for piecewise constant profiles
%  gsw_geo_strf_isopycnal                               - approximate isopycnal geostrophic streamfunction
%  gsw_geo_strf_isopycnal_pc                            - approximate isopycnal geostrophic streamfunction for
%                                                         piecewise constant profiles
%  gsw_geo_strf_Montgomery                              - Montgomery geostrophic streamfunction
%  gsw_geo_strf_Cunningham                              - Cunningham geostrophic streamfunction
%  gsw_geo_strf_steric_height                           - dynamic height anomaly divided by 9.7963 m s^-2
%  gsw_geo_strf_PISH                                    - pressure intergrated steric height
%  gsw_travel_time                                      - acoustic travel time
%  gsw_geostrophic_velocity                             - geostrophic velocity
%
% neutral versus isopycnal slopes and ratios
%  gsw_isopycnal_slope_ratio                            - ratio of the slopes of isopycnals on the SA-CT diagram 
%                                                         for p & p_ref
%  gsw_isopycnal_vs_ntp_CT_ratio                        - ratio of the gradient of Conservative Temperature
%                                                         in a potential density surface to that in the neutral 
%                                                         tangent plane
%  gsw_ntp_pt_vs_CT_ratio                               - ratio of gradients of potential temperature &
%                                                         Conservative Temperature in a neutral tangent plane
%                                                         (i.e. in a locally-referenced potential density surface)
%
% derivatives of entropy, CT and pt 
%  gsw_CT_first_derivatives                             - first derivatives of Conservative Temperature
%  gsw_CT_second_derivatives                            - second derivatives of Conservative Temperature
%  gsw_entropy_first_derivatives                        - first derivatives of entropy
%  gsw_entropy_second_derivatives                       - second derivatives of entropy
%  gsw_pt_first_derivatives                             - first derivatives of potential temperature
%  gsw_pt_second_derivatives                            - second derivatives of potential temperature
%
% seawater and ice properties at the freezing temperature
%  gsw_CT_freezing_poly                                 - Conservative Temperature freezing temperature of seawater (polynomial) 
%  gsw_t_freezing                                       - in-situ freezing temperature of seawater
%  gsw_t_freezing_poly                                  - in-situ freezing temperature of seawater (polynomial)
%  gsw_pot_enthalpy_ice_freezing                        - potential enthalpy of ice at which seawater freezes
%  gsw_pot_enthalpy_ice_freezing_poly                   - potential enthalpy of ice at which seawater freezes (polynomial)
%  gsw_SA_freezing_from_CT                              - Absolute Salinity of seawater at the freezing point (for given CT)
%  gsw_SA_freezing_from_CT_poly                         - Absolute Salinity of seawater at the freezing point (for given CT) (polynomial)
%  gsw_SA_freezing_from_t                               - Absolute Salinity of seawater at the freezing point (for given t)
%  gsw_SA_freezing_from_t_poly                          - Absolute Salinity of seawater at the freezing point (for given t) (polynomial)
%  gsw_pressure_freezing_CT                             - pressure of seawater at the freezing temperature (for given CT)
%  gsw_CT_freezing_first_derivatives                    - first derivatives of Conservative Temperature freezing temperature of seawater
%  gsw_CT_freezing_first_derivatives_poly               - first derivatives of Conservative Temperature freezing temperature of seawater (polynomial)
%  gsw_t_freezing_first_derivatives	first               - derivatives of in-situ freezing temperature of seawater
%  gsw_t_freezing_first_derivatives_poly                - first derivatives of in-situ freezing temperature of seawater (polynomial)
%  gsw_pot_enthalpy_ice_freezing_first_derivatives      - first derivatives o fpotential enthalpy of ice at which seawater freezes
%  gsw_pot_enthalpy_ice_freezing_first_derivatives_poly - first derivatives of potential enthalpy of ice at which seawater freezes (polynomial)
%  gsw_t_freezing_first_derivatives                     - first derivatives of in-situ freezing temperature of seawater
%  gsw_latentheat_melting                               - latent heat of melting of ice into seawater (isobaric melting enthalpy)
%
% thermodynamic interaction between ice and seawater
%  gsw_melting_ice_SA_CT_ratio                         - SA to CT ratio when ice melts in seawater
%  gsw_melting_ice_SA_CT_ratio_poly                    - SA to CT ratio when ice melts in seawater
%  gsw_melting_ice_equilibrium_SA_CT_ratio             - SA to CT ratio when ice melts into seawater near equilibrium
%  gsw_melting_ice_equilibrium_SA_CT_ratio_poly        - SA to CT ratio when ice melts into seawater near equilibrium
%  gsw_ice_fraction_to_freeze_seawater                 - ice mass fraction to freeze seawater
%  gsw_melting_ice_into_seawater                       - SA and CT when ice melts in seawater
%  gsw_frazil_ratios_adiabatic                         - ratios of SA, CT and P changes during frazil ice formation
%  gsw_frazil_ratios_adiabatic_poly                    - ratios of SA, CT and P changes during frazil ice formation
%  gsw_frazil_properties                               - SA, CT & ice mass fraction from bulk SA & bulk enthalpy
%  gsw_frazil_properties_potential                     - SA, CT & ice fraction from bulk SA & bulk potential enthalpy
%  gsw_frazil_properties_potential_poly                - SA, CT & ice fraction from bulk SA & bulk potential enthalpy (poly)
%
% thermodynamic interaction between sea ice and seawater
%  gsw_melting_seaice_SA_CT_ratio                      - SA to CT ratio when sea ice melts in seawater
%  gsw_melting_seaice_SA_CT_ratio                      - SA to CT ratio when sea ice melts in seawater (poly)
%  gsw_melting_seaice_equilibrium_SA_CT_ratio          - SA to CT ratio when sea ice melts into seawater equilibrium
%  gsw_melting_seaice_equilibrium_SA_CT_ratio_poly     - SA to CT ratio when sea ice melts into seawater equilibrium (poly)
%  gsw_melting_seaice_into_seawater                    - SA and CT when sea ice melts in seawater
%  gsw_seaice_fraction_to_freeze_seawater              - sea ice mass fraction to freeze seawater
%
% thermodynamic properties of ice Ih
%  gsw_rho_ice                                         - in-situ density of ice
%  gsw_alpha_wrt_t_ice                                 - thermal expansion coefficient of ice with respect to in-situ temperature
%  gsw_specvol_ice                                     - specific volume of ice
%  gsw_pressure_coefficient_ice                        - pressure coefficient of ice
%  gsw_sound_speed_ice                                 - sound speed of ice (compression waves)
%  gsw_kappa_ice                                       - isentropic compressibility of ice
%  gsw_kappa_const_t_ice                               - isothermal compressibility of ice
%  gsw_internal_energy_ice                             - internal energy of ice
%  gsw_enthalpy_ice                                    - enthalpy of ice
%  gsw_entropy_ice                                     - entropy of ice
%  gsw_cp_ice                                          - isobaric heat capacity of ice
%  gsw_chem_potential_water_ice                        - chemical potential of water in ice
%  gsw_Helmholtz_energy_ice                            - Helmholtz energy of ice
%  gsw_adiabatic_lapse_rate_ice                        - adiabatic lapse rate of ice
%  gsw_pt0_from_t_ice                                  - potential temperature of ice with reference pressure of 0 dbar
%  gsw_pt_from_t_ice                                   - potential temperature of ice
%  gsw_t_from_pt0_ice                                  - in-situ temperature from potential temperature of ice with p_ref of 0 dbar
%  gsw_t_from_rho_ice                                  - in-situ temperature from density of ice
%  gsw_pot_enthalpy_from_pt_ice                        - potential enthalpy from potential temperature of ice
%  gsw_pt_from_pot_enthalpy_ice                        - potential temperature from potential enthalpy of ice
%  gsw_pot_enthalpy_from_pt_ice_poly                   - potential enthalpy from potential temperature of ice (polynomial) 
%  gsw_pt_from_pot_enthalpy_ice_poly                   - potential temperature from potential enthalpy of ice (polynomial)
%  gsw_pot_enthalpy_from_specvol_ice                   - potential enthalpy from specific volume of ice 
%  gsw_specvol_from_pot_enthalpy_ice                   - specific volume from potential enthalpy of ice
%  gsw_pot_enthalpy_from_specvol_ice_poly              - potential enthalpy from specific volume of ice (polynomial) 
%  gsw_specvol_from_pot_enthalpy_ice_poly              - specific volume from potential enthalpy of ice (polynomial)
%
% isobaric evaporation enthalpy
%  gsw_latentheat_evap_CT                              - latent heat of evaporation of water from seawater (isobaric
%                                                        evaporation enthalpy) with CT as input temperature
%  gsw_latentheat_evap_t                               - latent heat of evaporation of water from seawater (isobaric
%                                                        evaporation enthalpy) with in-situ temperature, t, as input
%
% spiciness
%  gsw_spiciness0                                      - spiciness with reference pressure of 0 dbar
%  gsw_spiciness1                                      - spiciness with reference pressure of 1000 dbar
%  gsw_spiciness2                                      - spiciness with reference pressure of 2000 dbar
%  gsw_SA_CT_from_sigma0_spiciness0                    - SA & CT from given sigma and spiciness with p_ref of 0 dbar
%  gsw_SA_CT_from_sigma1_spiciness1                    - SA & CT from given sigma and spiciness with p_ref of 1000 dbar
%  gsw_SA_CT_from_sigma2_spiciness2                    - SA & CT from given sigma and spiciness with p_ref of 2000 dbar
%
% Planet Earth properties
%  gsw_f                                               - Coriolis parameter
%  gsw_grav                                            - gravitational acceleration
%  gsw_distance                                        - spherical earth distance between points in the ocean
%
% TEOS-10 constants
%  gsw_T0                                              - Celcius zero point; 273.15 K
%  gsw_P0                                              - one standard atmosphere; 101 325 Pa
%  gsw_SSO                                             - Standard Ocean Reference Salinity; 35.165 04 g/kg
%  gsw_uPS                                             - unit conversion factor for salinities; (35.165 04/35) g/kg
%  gsw_cp0                                             - the "specific heat" for use with CT; 3991.867 957 119 63 (J/kg)/K
%  gsw_C3515                                           - conductivity of SSW at SP=35, t_68=15, p=0; 42.9140 mS/cm
%  gsw_SonCl                                           - ratio of SP to Chlorinity; 1.80655 (g/kg)^-1
%  gsw_valence_factor                                  - valence factor of sea salt; 1.2452898
%  gsw_atomic_weight                                   - mole-weighted atomic weight of sea salt; 31.4038218... g/mol
%  
% dissolved gasses
%  gsw_Arsol                                           - argon solubility from SA and CT
%  gsw_Arsol_SP_pt                                     - argon solubility from SP and pt
%  gsw_Hesol                                           - helium solubility from SA and CT
%  gsw_Hesol_SP_pt                                     - helium solubility from SP and pt
%  gsw_Krsol                                           - krypton solubility from SA and CT
%  gsw_Krsol_SP_pt                                     - krypton solubility from SP and pt
%  gsw_N2sol                                           - nitrogen solubility from SA and CT
%  gsw_N2sol_SP_pt                                     - nitrogen solubility from SP and pt
%  gsw_Nesol                                           - neon solubility from SA and CT
%  gsw_Nesol_SP_pt                                     - neon solubility from SP and pt
%  gsw_O2sol                                           - oxygen solubility from SA and CT
%  gsw_O2sol_SP_pt                                     - oxygen solubility from SP and pt
%
% specific volume, density and enthalpy in terms of CT, based on the exact Gibbs function
%  gsw_specvol_CT_exact                                 - specific volume
%  gsw_alpha_CT_exact                                   - thermal expansion coefficient with respect to CT
%  gsw_beta_CT_exact                                    - saline contraction coefficientat constant CT
%  gsw_alpha_on_beta_CT_exact                           - alpha divided by beta
%  gsw_specvol_alpha_beta_CT_exact                      - specific volume, thermal expansion & saline contraction coefficients
%  gsw_specvol_first_derivatives_CT_exact               - first derivatives of rho 
%  gsw_specvol_second_derivatives_CT_exact              - second derivatives of specific volume
%  gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact  - first derivatives of specific volume with respect to enthalpy
%  gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact - second derivatives of specific volume with respect to enthalpy
%  gsw_specvol_anom_CT_exact                            - specific volume anomaly
%  gsw_specvol_anom_standard_CT_exact                   - specific volume anomaly relative to SSO & 0C
%  gsw_specvol_diff_CT_exact                            - specific volume difference between two pressures
%  gsw_rho_CT_exact                                     - in-situ density and potential density
%  gsw_rho_alpha_beta_CT_exact                          - in-situ density, thermal expansion & saline contraction coefficient
%  gsw_rho_first_derivatives_CT_exact                   - first derivatives of density
%  gsw_rho_second_derivatives_CT_exact                  - second derivatives of density
%  gsw_rho_first_derivatives_wrt_enthalpy_CT_exact      - first derivatives of density with respect to enthalpy
%  gsw_rho_second_derivatives_wrt_enthalpy_CT_exact     - second derivatives of density with respect to enthalpy
%  gsw_sigma0_CT_exact                                  - sigma0 with reference pressure of 0 dbar
%  gsw_sigma1_CT_exact                                  - sigma1 with reference pressure of 1000 dbar
%  gsw_sigma2_CT_exact                                  - sigma2 with reference pressure of 2000 dbar
%  gsw_sigma3_CT_exact                                  - sigma3 with reference pressure of 3000 dbar
%  gsw_sigma4_CT_exact                                  - sigma4 with reference pressure of 4000 dbar
%  gsw_cabbeling_CT_exact                               - cabbeling coefficient
%  gsw_thermobaric_CT_exact                             - thermobaric coefficient
%  gsw_enthalpy_CT_exact                                - enthalpy
%  gsw_enthalpy_diff_CT_exact                           - difference of enthalpy between two pressures
%  gsw_dynamic_enthalpy_CT_exact                        - dynamic enthalpy
%  gsw_enthalpy_first_derivatives_CT_exact              - first derivatives of enthalpy
%  gsw_enthalpy_second_derivatives_CT_exact             - second derivatives of enthalpy
%  gsw_sound_speed_CT_exact                             - sound speed
%  gsw_kappa_CT_exact                                   - isentropic compressibility
%  gsw_internal_energy_CT_exact                         - internal energy
%  gsw_internal_energy_first_derivatives_CT_exact       - first derivatives of internal energy
%  gsw_internal_energy_second_derivatives_CT_exact      - second derivatives of internal energy
%  gsw_CT_from_enthalpy_exact                           - Conservative Temperature from enthalpy
%  gsw_SA_from_rho_CT_exact                             - Absolute Salinity from density
%  gsw_CT_from_rho_exact                                - Conservative Temperature from density
%  gsw_CT_maxdensity_exact                              - Conservative Temperature of maximum density of seawater
%
% Laboratory functions, for use with densimeter measurements
%  gsw_SA_from_rho_t_exact                              - Absolute Salinity from density
%  gsw_deltaSA_from_rho_t_exact                         - Absolute Salinity Anomaly from density
%  gsw_rho_t_exact                                      - in-situ density
%
% basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
%  gsw_specvol_t_exact                                  - specific volume
%  gsw_alpha_wrt_CT_t_exact                             - thermal expansion coefficient with respect to 
%                                                         Conservative Temperature.
%  gsw_alpha_wrt_pt_t_exact                             - thermal expansion coefficient with respect to 
%                                                         potential temperature
%  gsw_alpha_wrt_t_exact                                - thermal expansion coefficient with respect to 
%                                                         in-situ temperature
%  gsw_beta_const_CT_t_exact                            - saline contraction coefficient at constant 
%                                                         Conservative Temperature.
%  gsw_beta_const_pt_t_exact                            - saline contraction coefficient at constant 
%                                                         potential temperature
%  gsw_beta_const_t_exact                               - saline contraction coefficient at constant 
%                                                         in-situ temperature
%  gsw_specvol_anom_t_exact                             - specific volume anomaly
%  gsw_rho_t_exact                                      - in-situ density
%  gsw_pot_rho_t_exact                                  - potential density
%  gsw_sigma0_pt0_exact                                 - sigma0 from pt0 with reference pressure of 0 dbar
%  gsw_enthalpy_t_exact                                 - enthalpy
%  gsw_dynamic_enthalpy_t_exact                         - dynamic enthalpy
%  gsw_CT_first_derivatives_wrt_t_exact	first           - derivatives of Conservative Temperature with respect to t
%  gsw_enthalpy_first_derivatives_wrt_t_exact           - first derivatives of enthalpy with respect to t
%  gsw_sound_speed_t_exact                              - sound speed
%  gsw_kappa_t_exact                                    - isentropic compressibility
%  gsw_kappa_const_t_exact                              - isothermal compressibility
%  gsw_internal_energy_t_exact                          - internal energy
%  gsw_SA_from_rho_t_exact                              - Absolute Salinity from density
%  gsw_t_from_rho_exact                                 - in-situ temperature from density
%  gsw_t_maxdensity_exact                               - in-situ temperature of maximum density of seawater
%  gsw_cp_t_exact                                       - isobaric heat capacity
%  gsw_isochoric_heat_cap_t_exact                       - isochoric heat capacity
%  gsw_chem_potential_relative_t_exact                  - relative chemical potential
%  gsw_chem_potential_water_t_exact                     - chemical potential of water in seawater
%  gsw_chem_potential_salt_t_exact                      - chemical potential of salt in seawater
%  gsw_Gibbs_energy_t_exact                             - Gibbs energy
%  gsw_Helmholtz_energy_t_exact                         - Helmholtz energy
%  gsw_osmotic_coefficient_t_exact                      - osmotic coefficient of seawater
%  gsw_osmotic_pressure_t_exact                         - osmotic pressure of seawater
%
% Library functions of the GSW toolbox (internal functions; not intended to be called by users) 
%  (The GSW functions above call the following library functions.)
%  gsw_gibbs                                           - the TEOS-10 Gibbs function and its derivatives
%  gsw_gibbs_ice                                       - the TEOS-10 Gibbs function of ice and its derivatives
%  gsw_SAAR                                            - Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
%  gsw_Fdelta                                          - ratio of Absolute to Preformed Salinity, minus 1
%  gsw_deltaSA_atlas                                   - Absolute Salinity Anomaly atlas value (excluding the Baltic Sea)
%  gsw_SA_from_SP_Baltic                               - calculates Absolute Salinity in the Baltic Sea
%  gsw_SP_from_SA_Baltic                               - calculates Practical Salinity in the Baltic Sea
%  gsw_infunnel                                        - "oceanographic funnel" check for the 75-term equation
%  gsw_entropy_part                                    - entropy minus the terms that are a function of only SA
%  gsw_entropy_part_zerop                              - entropy_part evaluated at 0 dbar
%  gsw_quadprog                                        - quadratic solver for water column stabilisation
%  gsw_wiggliness                                      - amount of variation in a cast
%  gsw_data_interp                                     - Barker & McDougall (2020) MR-pchip interpolation
%  gsw_interp_ref_cast                                 - linearly interpolates the reference cast
%  gsw_linear_interp_SA_CT                             - linearly interpolates (SA,CT,p) to the desired p
%  gsw_pchip_interp_SA_CT                              - pchip interpolation of (SA,CT,p) to the desired p
%  gsw_rr68_interp_SA_CT                               - Reininger & Ross (1968) interpolation of (SA,CT,p) to the desired p
%  gsw_spline_interp_SA_CT                             - spline interpolation of (SA,CT,p) to the desired p
%  gsw_wiggliness                                      - amount of variation in a cast
%  gsw_gibbs_pt0_pt0                                   - gibbs(0,2,0,SA,t,0)
%  gsw_gibbs_pt0_pt0                                   - gibbs(0,2,0,SA,t,0)
%  gsw_gibbs_ice_part_t                                - part of gibbs_ice(1,0,t,p)
%  gsw_gibbs_ice_pt0                                   - part of gibbs_ice(1,0,pt0,0)
%  gsw_specvol_SSO_0                                   - specvol(35.16504,0,p)
%  gsw_enthalpy_SSO_0                                  - enthalpy(35.16504,0,p)
%  gsw_Hill_ratio_at_SP2                               - Hill ratio at a Practical Salinity of 2
% 
%  The GSW data set.
%  gsw_data_v3_0        - This file contains:
%                          (1) the global data set of Absolute Salinity Anomaly Ratio,
%                          (2) the global data set of Absolute Salinity Anomaly Ref.,                                    
%                          (3) a reference cast (for the isopycnal streamfunction), 
%                          (4) two reference casts that are used by gsw_demo 
%                          (5) three vertical profiles of (SP, t, p) at known long & lat, plus
%                              the outputs of all the GSW functions for these 3 profiles, and
%                              the required accuracy of all these outputs.
%