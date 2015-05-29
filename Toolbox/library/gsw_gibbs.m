% gsw_gibbs                                Gibbs energy and its derivatives
%==========================================================================
%
% USAGE:
%   gibbs = gsw_gibbs(ns,nt,np,SA,t,p)
%
% DESCRIPTION:
%  Calculates specific Gibbs energy and its derivatives up to order 3 for
%  seawater.  The Gibbs function for seawater is that of TEOS-10
%  (IOC et al., 2010), being the sum of IAPWS-08 for the saline part and
%  IAPWS-09 for the pure water part.  These IAPWS releases are the
%  officially blessed IAPWS descriptions of Feistel (2008) and the pure
%  water part of Feistel (2003).  Absolute Salinity, SA, in all of the GSW
%  routines is expressed on the Reference-Composition Salinity Scale of
%  2008 (RCSS-08) of Millero et al. (2008).
%
% INPUT:
%  ns  =  order of SA derivative                     [ integers 0, 1 or 2 ]
%  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
%  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%
%  SA, t and p need to have the same dimensions.
%
% OUTPUT:
%  gibbs  =  Specific Gibbs energy or its derivatives.
%
%            The Gibbs energy (when ns = nt = np = 0) has units of:
%                                                                  [ J/kg ]
%            The Absolute Salinity derivatives are output in units of:
%                                                   [ (J/kg) (g/kg)^(-ns) ]
%            The temperature derivatives are output in units of:
%                                                      [ (J/kg) (K)^(-nt) ]
%            The pressure derivatives are output in units of:
%                                                     [ (J/kg) (Pa)^(-np) ]
%            The mixed derivatives are output in units of:
%                              [ (J/kg) (g/kg)^(-ns) (K)^(-nt) (Pa)^(-np) ]
%  Note. The derivatives are taken with respect to pressure in Pa, not
%    withstanding that the pressure input into this routine is in dbar.
%
% AUTHOR:
%  David Jackett, Paul Barker and Trevor McDougall     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  Feistel, R., 2003: A new extended Gibbs thermodynamic potential of
%   seawater,  Progr. Oceanogr., 58, 43-114.
%
%  Feistel, R., 2008: A Gibbs function for seawater thermodynamics
%   for -6 to 80°C and salinity up to 120 g kg–1, Deep-Sea Res. I,
%   55, 1639-1671.
%
%  IAPWS, 2008: Release on the IAPWS Formulation 2008 for the
%   Thermodynamic Properties of Seawater. The International Association
%   for the Properties of Water and Steam. Berlin, Germany, September
%   2008, available from http://www.iapws.org.  This Release is referred
%   to as IAPWS-08.
%
%  IAPWS, 2009: Supplementary Release on a Computationally Efficient
%   Thermodynamic Formulation for Liquid Water for Oceanographic Use.
%   The International Association for the Properties of Water and Steam.
%   Doorwerth, The Netherlands, September 2009, available from
%   http://www.iapws.org.  This Release is referred to as IAPWS-09.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.6 and appendices A.6,  G and H of this TEOS-10 Manual.
%
%  Millero, F.J., R. Feistel, D.G. Wright, and T.J. McDougall, 2008:
%   The composition of Standard Seawater and the definition of the
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================