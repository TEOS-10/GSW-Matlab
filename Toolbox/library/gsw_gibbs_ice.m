% gsw_gibbs_ice                     Gibbs energy of ice and its derivatives
% =========================================================================
%
% USAGE:
%  gibbs_ice = gsw_gibbs_ice(nt,np,t,p)
%
% DESCRIPTION:
%  Ice specific Gibbs energy and derivatives up to order 2.
%
% INPUT:
%  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
%  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%   
% OUTPUT:
%  gibbs_ice = Specific Gibbs energy of ice or its derivatives.
%            The Gibbs energy (when nt = np = 0) has units of:     [ J/kg ]
%            The temperature derivatives are output in units of: 
%                                                      [ (J/kg) (K)^(-nt) ]
%            The pressure derivatives are output in units of:
%                                                     [ (J/kg) (Pa)^(-np) ]
%            The mixed derivatives are output in units of:
%                                           [ (J/kg) (K)^(-nt) (Pa)^(-np) ]
%  Note. The derivatives are taken with respect to pressure in Pa, not
%    withstanding that the pressure input into this routine is in dbar.
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IAPWS, 2009: Revised release on the Equation of State 2006 for H2O Ice 
%   Ih. The International Association for the Properties of Water and 
%   Steam. Doorwerth, The Netherlands, September 2009.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See appendix I.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================