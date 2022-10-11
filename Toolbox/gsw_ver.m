function gsw_ver

% gsw_ver                          GSW Oceanographic Toolbox version number
%==========================================================================
%
% USAGE:  
%  gsw_ver
%
% DESCRIPTION:
%  This function displays the version number of the GSW Oceanographic
%  Toolbox.
%
% AUTHOR:  
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) and section 3.27 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

gsw_data = 'gsw_data_v3_0.mat';
gsw_data_file = which(gsw_data);
load(gsw_data_file,'version_number','version_date');

disp('====================================================================')
disp(' ')
disp('  Gibbs SeaWater (GSW) Oceanographic Toolbox')
disp(['  Version ',version_number])
disp('  This toolbox uses the Absolute Salinity Anomaly Ratio ')
disp(['  look-up table dataset ',gsw_data(1:(end-4)),' (',version_date,')'])
disp(' ')
disp('====================================================================')

end
