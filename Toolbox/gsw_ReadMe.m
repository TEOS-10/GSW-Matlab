function gsw_ReadMe

%%
% ***************************************************************
% This is the gsw (Gibbs sea water) library implemented in Matlab
% ***************************************************************
%
%   It consists of two parts
%
%   (i) the actual Gibbs function and derived major oceanographic variables.

%   The Gibbs function is the sum of pure water and saline components defined by
%    
%       Feistel, R., 2003: A new extended Gibbs thermodynamic potential of seawater,  
%               Progr. Oceanogr., 58, 43-114
%
%   for pure water, and for the saline part
%
%       Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80°C
%               and salinity up to 120 g/kg, Deep-Sea Res. I, 55, 1639-1671, and  
%     
%       IAPWS 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic Properties
%               of Seawater, The International Association for the Properties of Water
%               and Steam, Berlin, Germany, September 2008,
%               available at http://www.iapws.org.
%  	
%   (ii) the conversion of Practical Salinity to Absolute Salinity as described in 
%
%       McDougall, T.J., Jackett. D.R. and Millero, F.J., 2009: An algorithm for estimating 
%               Absolute Salinity in the global ocean, Ocean Science, submitted.  
%     
%  	
% Implementation by David Jackett
% January 2009
