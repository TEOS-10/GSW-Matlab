function data = gsw_resize(data,sz_SA)

% gsw_resize                                                   resizes data
%==========================================================================
%
% USAGE:  
%  data = gsw_resize(data,sz_SA)
%
% DESCRIPTION:
%  This function resizes the input data such that it has the same size as
%  salinity.  
%
% INPUT:
%  data   =  input data
%  sz_SA  =  size of Absolute Salinity                                        
%
%  data is resized if it has up to 3 dimensions.
%
% OUTPUT:
%  data  =  data resized such that it has the same size as Absolute 
%           Salinity sz_SA
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
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

sz_data = size(data);

dims_SA = length(sz_SA);
dims_data = length(sz_data);

if dims_SA == dims_data & sum(abs(sz_SA - sz_data)) == 0
    if sum(abs(sz_SA - sz_data)) == 0
        % ok
    end
    
else  % resizing is needed
    if dims_SA == 2
        if sz_data(1) == 1 & sz_data(2) == 1
            data = data*ones(sz_SA(1),sz_SA(2));
        elseif sz_SA(1) == sz_data(1) & sz_data(2) == 1
            data = data(ones(1,sz_SA(2)), :);
        elseif sz_data(1) == 1 & sz_SA(2) == sz_data(2)
            data = data(:,ones(1,sz_SA(1)));
        elseif sz_data(1) == 1 & sz_SA(1) == sz_data(2)
            data = data.';
            data = data(ones(1,sz_SA(2)),:);       
        elseif (sz_SA(2) == sz_data(1)) & (sz_data(2) == 1)
            data = data.';
            data = data(:,ones(1,sz_SA(1)));       
        elseif sz_SA(1) == sz_data(2) & sz_SA(2) == sz_data(1) & sz_SA(1) ~= sz_SA(2)
            data = data.';
        elseif sz_SA(1) == sz_data(1) & sz_SA(2) == sz_data(2)
            % ok
        else
            error('gsw_resize: Inputs array dimensions arguments do not agree')
        end
               
    elseif dims_SA == 3
        sz_data = size(squeeze(data));
        if length(sz_data) == 2

            % scalar
            if sz_data(1) == 1 & sz_data(2) == 1
                data = data*ones(sz_SA(1),sz_SA(2),sz_SA(3));
                
            % vector
            elseif (sz_data(1) == 1 | sz_data(2) == 1)  & sz_data(1) ~= sz_data(2)
                data = data(:);
                lp = size(data);
                
                if sz_SA(1) == lp(1) & sz_SA(2) ~= lp(1) & sz_SA(3) ~= lp(1)
                    data_dummy(:,1,1) = data;
                    data = repmat(data_dummy,[1,sz_SA(2),sz_SA(3)]);
                elseif sz_SA(1) ~= lp(1) & sz_SA(2) == lp(1) & sz_SA(3) ~= lp(1)
                    data_dummy(1,:,1) = data;
                    data = repmat(data_dummy,[sz_SA(1),1,sz_SA(3)]);
                elseif sz_SA(1) ~= lp(1) & sz_SA(2) ~= lp(1) & sz_SA(3) == lp(1)
                    data_dummy(1,1,:) = data;
                    data = repmat(data_dummy,[sz_SA(1),sz_SA(2),1]);
                else
                    error('gsw_resize: Inputs array dimensions arguments do not agree')
                end
                
            % matrix
            elseif  sz_data(1) == sz_SA(1) & sz_data(1) ~= sz_SA(2) & sz_data(1) ~= sz_SA(3) & sz_data(2) ~= sz_SA(1) & sz_data(2) == sz_SA(2) & sz_data(2) ~= sz_SA(3)
                data_dummy(:,:,1) = squeeze(data);
                data = repmat(data_dummy,[1,1,sz_SA(3)]);
            elseif  sz_data(1) == sz_SA(1) & sz_data(1) ~= sz_SA(2) & sz_data(1) ~= sz_SA(3) & sz_data(2) ~= sz_SA(1) & sz_data(2) ~= sz_SA(2) & sz_data(2) == sz_SA(3)
                data_dummy(:,1,:) = squeeze(data);
                data = repmat(data_dummy,[1,sz_SA(2),1]);
            elseif  sz_data(1) ~= sz_SA(1) & sz_data(1) == sz_SA(2) & sz_data(1) ~= sz_SA(3) & sz_data(2) == sz_SA(1) & sz_data(2) ~= sz_SA(2) & sz_data(2) ~= sz_SA(3)
                data = (squeeze(data)).';
                data_dummy(:,:,1) = squeeze(data);
                data = repmat(data_dummy,[1,1,sz_SA(3)]);
            elseif  sz_data(1) ~= sz_SA(1) & sz_data(1) == sz_SA(2) & sz_data(1) ~= sz_SA(3) & sz_data(2) ~= sz_SA(1) & sz_data(2) ~= sz_SA(2) & sz_data(2) == sz_SA(3)
                data_dummy(1,:,:) = squeeze(data);
                data = repmat(data_dummy,[sz_SA(1),1,1]);
            elseif  sz_data(1) ~= sz_SA(1) & sz_data(1) ~= sz_SA(2) & sz_data(1) == sz_SA(3) & sz_data(2) == sz_SA(1) & sz_data(2) ~= sz_SA(2) & sz_data(2) ~= sz_SA(3)
                data = (squeeze(data)).';
                data_dummy(:,1,:) = data;
                data = repmat(data_dummy,[1,sz_SA(2),1]);
            elseif  sz_data(1) ~= sz_SA(1) & sz_data(1) ~= sz_SA(2) & sz_data(1) == sz_SA(3) & sz_data(2) ~= sz_SA(1) & sz_data(2) == sz_SA(2) & sz_data(2) ~= sz_SA(3)
                data = (squeeze(data)).';
                data_dummy(1,:,:) = data;
                data = repmat(data_dummy,[sz_SA(1),1,1]);
            else
                error('gsw_resize: Inputs array dimensions arguments do not agree')
            end
            
        elseif length(sz_data) == 3
            if sum(abs(sz_SA - sz_data)) ~= 0
                error('gsw_resize: Inputs array dimensions arguments do not agree')
            else
                % ok
            end
        else
            error('gsw_resize: Inputs array dimensions arguments do not agree')
        end
        
    else % there are more than 3 dimensions, the dimensions must be equal
        if sum(abs(sz_SA - sz_data)) ~= 0
            error('gsw_resize: Inputs array dimensions arguments do not agree')
        else
            % ok
        end
    end
end

end
