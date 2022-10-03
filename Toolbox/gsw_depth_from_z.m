function depth = gsw_depth_from_z(z)

% gsw_depth_from_z                                     depth from height, z
%==========================================================================
%
% USAGE:  
%  depth = gsw_depth_from_z(z)
%
% DESCRIPTION:
%  Calculates depth from height, z.  Note that in general height is
%  negative in the ocean.   
%
% INPUT:
%  z  =  height                                                       [ m ]
%
% OUTPUT:
%  depth  =  depth                                                    [ m ]
%
% AUTHOR:  
%  Winston                                              [ god@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%   
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_depth_from_z: Requires one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

depth = -z;

end
