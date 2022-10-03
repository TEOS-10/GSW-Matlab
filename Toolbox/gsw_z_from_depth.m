function z = gsw_z_from_depth(depth)

% gsw_z_from_depth                                    height, z, from depth
%==========================================================================
%
% USAGE:  
%  z = gsw_z_from_depth(depth)
%
% DESCRIPTION:
%  Calculates height, z, from depth.  Note that in general height is
%  negative in the ocean.   
%
% INPUT:
%  depth  =  depth                                                    [ m ]
%
% OUTPUT:
%  z  =  height                                                       [ m ]
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
   error('gsw_z_from_depth: Requires one input')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

z = -depth;

end
