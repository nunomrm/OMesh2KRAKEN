function [del,obj] = setProj(obj,proj,projtype,insert,bbox)
    % [del,obj] = setProj(obj,proj,projtype,insert)
    % kjr generic function to parse projected space options
    % returns del flag to delete overlapping elements when plotting
    % global models.
    % if insert = 1, then automatically insert the global m_proj variables
    % into the msh obj. insert is 0 by default. 

    if nargin < 4
        insert = 0;
    end
    
    if nargin < 5
        bbox = [];
    end
    
    % process bounds of mesh (or supply your own)
    if isempty(bbox)
        lon_mi = min(obj.p(:,1)); lon_ma = max(obj.p(:,1));
        lat_mi = min(obj.p(:,2)); lat_ma = max(obj.p(:,2));
    else
        lon_mi = bbox(1,1); lon_ma = bbox(1,2);
        lat_mi = bbox(2,1); lat_ma = bbox(2,2);
    end
    lat_mea = mean(obj.p(:,2)); lon_mea = mean(obj.p(:,1));
    
    % some defaults
    rad = 100; rot = 15;
    del = 0 ;

    % process the projtype as a varargin type
    I = find(strcmp(projtype,'long'));
    if ~isempty(I)
        if length(projtype{I+1}) == 1
            lon_mea = projtype{I+1};
        else
            lon_mi = projtype{I+1}(1); lon_ma = projtype{I+1}(2);
        end
    end
    I = find(strcmp(projtype,'lat'));
    if ~isempty(I)
        if length(projtype{I+1}) == 1
            lat_mea = projtype{I+1};
        else
            lat_mi = projtype{I+1}(1); lat_ma = projtype{I+1}(2);
        end
    end
    I = find(strcmp(projtype,'rad'));
    if ~isempty(I)
        rad = projtype{I+1};
    end
    I = find(strcmp(projtype,'rot'));
    if ~isempty(I)
        rot = projtype{I+1};
    end
    if proj == 0
        % normal geographic coordinates
        m_proj('equi','lat',[lat_mi lat_ma],'long',[lon_mi lon_ma]) ;
        del = 1;
    else
        if ~ischar(projtype)
            projtype = projtype{1};
        end
        projtype = lower(projtype);
        if ~isempty(regexp(projtype,'ste'))
            % Special treatment of Stereographic projection
            if lat_ma < 0
                % center Antarctica
                m_proj(projtype,'lat',-90,...
                      'long',0.5*(lon_mi+lon_ma),...
                      'radius',min(lat_ma+90,180),'rot',rot);
            else
                % center Arctic
                m_proj(projtype,'lat',90,...
                      'long',0.5*(lon_mi+lon_ma),...
                      'radius',min(90-lat_mi,180),'rot',rot);
            end
            m_proj('get') ;
        elseif  ~isempty(regexp(projtype,'ort')) || ...
                ~isempty(regexp(projtype,'gno')) || ...
                ~isempty(regexp(projtype,'azi')) || ...
                ~isempty(regexp(projtype,'sat'))
            m_proj(projtype,'lat',lat_mea,'long',lon_mea,...
                   'radius',rad,'rot',rot);
            m_proj('get') ;
            del = 1;
        elseif ~isempty(regexp(projtype,'obl')) 
            % Oblique Mercator projection
            asp = (lon_ma-lon_mi)/(lat_ma - lat_mi);
            dir = 'hor';
            if asp > 1 
                asp = 1./asp; dir = 'ver';
            end
            m_proj(projtype,'lon',[lon_mi lon_ma],...
                            'lat',[lat_mi lat_ma],...
                            'aspect',asp,'dir',dir) ;
            m_proj('get') ;
            del = 1;
        else
            % Cylindrical, Conic or Global type projections
            del = 1;
            m_proj(projtype,'lon',[lon_mi lon_ma],...
                            'lat',[lat_mi lat_ma]) ;
            m_proj('get') ;
        end
    end
    if insert
        global MAP_PROJECTION MAP_COORDS MAP_VAR_LIST
        obj.proj   = MAP_PROJECTION ; 
        obj.coord  = MAP_COORDS ; 
        obj.mapvar = MAP_VAR_LIST ; 
    end
end
