clear all, close all, clc

SP = 35.52764437773385
SA = 35.7
t = 25.5
p = 1023
pr = 0,
longs = 201
lats = -21,

disp('****************************************************************')    
first_column = '                   computed value' 
second_column = '                   expected computed value'
disp('****************************************************************')       

ASal = [gsw_ASal(SP,p,longs,lats), 35.7]

PSal_from_ASal = [gsw_PSal_from_ASal(SA,p,longs,lats), 35.5276443777339]


alpha_t = [gsw_alpha_t(SA,t,p), 0.000309837839319264]

beta_t = [gsw_beta_t(SA,t,p), 0.000725729797838666]

cp = [gsw_cp(SA,t,p), 3974.42541259729]
           
ctmp = [gsw_ctmp(SA,t), 25.4805463842239] 
           
ptmp0_from_ctmp = [gsw_ptmp0_from_ctmp(SA,ctmp(1)), 25.5] 

density = [gsw_dens(SA,t,p),1027.95249315662]
           
enthalpy = [gsw_enthalpy(SA,t,p), 110776.712408975]
                      
entropy = [gsw_entropy(SA,t,p), 352.81879771528]
                      
kappa = [gsw_kappa(SA,t,p), 4.03386268546478e-6]
                                 
kappa_t = [gsw_kappa_t(SA,t,p), 4.10403794615135e-6]
                                 
pden = [gsw_pden(SA,t,p,pr), 1023.66254941185]       
                                 
ptmp = [gsw_ptmp(SA,t,p,pr), 25.2720983155409]       
                                 
ptmp_inverse = [gsw_ptmp(SA,ptmp(1),pr,p), 25.5] 

specvol = [gsw_specvol(SA,t,p), 0.0009728076021579713] 
           
svel = [gsw_svel(SA,t,p), 1552.93372863425] 

           