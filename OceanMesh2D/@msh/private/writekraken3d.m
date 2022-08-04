function outname = writekraken3d(outfiname,tri,VX,B,fname_env,TS_data,params_env,coords_or,params_flp,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                    %
% tri: triangle (nelem,3);                                  %
% VX(:,1): longitude (nnode,1);                             %
% VX(:,2): latitude (nnode,1);                              %
% B: depth (nnode,1);                                       %
% outfiname: mesh name                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_env=1; % flag to write ENV files (0: no; 1: yes)
flag_flp=1; % flag to write the FLP file (0: no; 1: yes)

if (nargin >= 10)
   flags=varargin{1};
   flag_env = flags(1);
   flag_flp = flags(2);
end

nnode=length(VX(:,1));
lon=VX(:,1);
lat=VX(:,2);

while flag_env==1
    
    disp('write ENV files for KRAKEN and a FLP file');
       
    % T-S data extraction
    lo_d=TS_data.lon;
    la_d=TS_data.lat;
    T=TS_data.T;
    S=TS_data.S;
    z_TS=TS_data.z;
    x1=min(lon); x2=max(lon); y1=min(lat); y2=max(lat);
    x_f=lo_d>=x1-0.1&lo_d<=x2+0.1;
    y_f=la_d>=y1-0.1&la_d<=y2+0.1;
    lo_d=lo_d(x_f);
    la_d=la_d(y_f);
    [LO_D,LA_D]=meshgrid(lo_d,la_d);
    LO_D=LO_D';
    LA_D=LA_D';
    
    % Calculation of SSP in the unstructured mesh from T-S data
    idxx=[];
    N=nnode;
    Nz=length(z_TS);
    Nx_d=length(lo_d);
    Ny_d=length(la_d);
    S_n=zeros(N,Nz);
    T_n=zeros(N,Nz);
    D_n=zeros(N,Nz);
    SSP_n=zeros(N,Nz);
    
    for i_z=1:Nz
        S_z=squeeze(S(x_f,y_f,i_z,1));
        T_z=squeeze(T(x_f,y_f,i_z,1));
        c=1;
        size(S_z)
        length(lo_d)
        length(la_d)
        size(LO_D)
        S_zz=zeros(Nx_d*Ny_d,1);
        T_zz=zeros(Nx_d*Ny_d,1);
        lo_zz=zeros(Nx_d*Ny_d,1);
        la_zz=zeros(Nx_d*Ny_d,1);
        for i=1:Nx_d
            for j=1:Ny_d
                S_zz(c)=S_z(i,j);
                T_zz(c)=T_z(i,j);
                lo_zz(c)=LO_D(i,j);
                la_zz(c)=LA_D(i,j);
                c=c+1;
            end
        end
        F=scatteredInterpolant(double(lo_zz),double(la_zz),S_zz);
        S_n(:,i_z) = F(lon,lat) ;
        F=scatteredInterpolant(double(lo_zz),double(la_zz),T_zz);
        T_n(:,i_z) = F(lon,lat) ;
    
        % sounspeed profiles
        SSP_n(:,i_z)=sndspd(S_n(:,i_z),T_n(:,i_z),z_TS(i_z));
        if sum(isnan(SSP_n(:,i_z)))~=N
            for i_m=1:N
                lon_search=lon;
                lat_search=lat;
                lon_search2=lon;
                lat_search2=lat;
                lon_search3=lon;
                lat_search3=lat;
                if B(i_m)<=z_TS(i_z) & z_TS(i_z)>5
                    SSP_n(i_m,i_z)=NaN;
                end
                if z_TS(i_z)<=B(i_m) & (isnan(SSP_n(i_m,i_z)) || SSP_n(i_m,i_z)<1500)
                    while 1
                        [lo,la,idx]=find_nearest_point(lon_search3,lat_search3,lon(i_m),lat(i_m));
                        if isnan(SSP_n(idx,i_z))
                            lon_search3(idx)=NaN;
                            lat_search3(idx)=NaN;
                        else
                            idxx(end+1)=idx;
                            SSP_n(i_m,i_z)=SSP_n(idx,i_z);
                            break
                        end
                    end
                end
                if isnan(SSP_n(i_m,1))
                    while 1
                        [lo,la,idx]=find_nearest_point(lon_search,lat_search,lon(i_m),lat(i_m));
                        if isnan(SSP_n(idx,1))
                            lon_search(idx)=NaN;
                            lat_search(idx)=NaN;
                        else
                            SSP_n(i_m,1)=SSP_n(idx,1);
                            if isnan(SSP_n(i_m,2))
                                SSP_n(i_m,2)=SSP_n(i_m,1);
                            end
                            break
                        end
                    end
                end
                if isnan(SSP_n(i_m,2))
                    while 1
                        [lo,la,idx]=find_nearest_point(lon_search2,lat_search2,lon(i_m),lat(i_m));
                        if isnan(SSP_n(idx,1))
                            lon_search2(idx)=NaN;
                            lat_search2(idx)=NaN;
                        else
                            SSP_n(i_m,2)=SSP_n(i_m,1);
                            break
                        end
                    end
                end
                
            end
        end
                
    end
    
    
    % remove NaN layers of SSP
    c=0;
    size(SSP_n)
    Nz
    for i_z=1:Nz
        if sum(isnan(SSP_n(:,i_z-c)))==N
            SSP_n(:,i_z-c)=[];
            z_TS(i_z-c)=[];
            c=c+1;
        end
    end
    
    D=B; % topo-bathymetry
    disp('Writing .env files')
    cInt.Low = params_env{end-3};
    cInt.High = params_env{end-2};
    model    = 'KRAKEN';
    Receiverfin = params_env{8}; 
    deltassp=params_env{end-1};
    sourcedepth=params_env{end-5};
    freqT    = params_env{1};
    SSP.NMedia = params_env{2};
    SSP.N      = [ params_env{6:7} ];
    SSP.sigma  = [ params_env{4:5} ]; %Hrms
    z_fin=1:deltassp:max(B);
    SSP_fin = zeros(N,length(z_fin)).*NaN;
    
    
    for idx=1:N 
        envfil   = sprintf('%s_%06d',fname_env,idx);
        TitleEnv = envfil;
        zi=0:deltassp:D(idx);
    
        %Water
        if D(idx)>=deltassp 
            if round(10.^2*zi(length(zi)))/10.^2==round(10.^2*D(idx))/10.^2
                ftest(idx)=1;   
            else
                ftest(idx)=0;
                zi(length(zi)+1)=D(idx);
            end
            SSP_x=interp1(z_TS',SSP_n(idx,:),zi','spline');
            SSP_fin(idx,1:length(SSP_x))=interp1(z_TS',SSP_n(idx,:),zi','spline');
        elseif D(idx)<deltassp
            if zi<=1e-2
                zi=[0 1e-2];
            else
                zi=[0 D(idx)];
            end
        end
        SSP.depth  = [0 zi(end) Receiverfin ];        
        WalphaR=interp1(z_TS',SSP_n(idx,:),zi','spline');
        
        f_z=find(zi>max(z_TS));
        if length(f_z)>0
            WalphaR(f_z)=ones(1,length(f_z)).*WalphaR(f_z(1)-1);
        end
    
        WbetaR=ones(size(zi,2),1)*0;
        Wrho=ones(size(zi,2))*1;
        WalphaI=ones(size(zi,2))*0;
        WbetaI=ones(size(zi,2))*0;      
        SSP.raw( 1 ).z      = zi;
        SSP.raw( 1 ).alphaR = WalphaR;
        SSP.raw( 1 ).betaR  = WbetaR;
        SSP.raw( 1 ).rho    = Wrho;
        SSP.raw( 1 ).alphaI = WalphaI;
        SSP.raw( 1 ).betaI  = WbetaI;                  
        %Bottom
        SSP.raw( 2 ).z      = [ zi(end) Receiverfin ];
        SSP.raw( 2 ).alphaR = [ params_env{end} params_env{end} ]; % speed
        Bdry.Bot.HS.alphaR = 0;
        SSP.raw( 2 ).betaR  = [ 0   0 ];
        Bdry.Bot.HS.betaR= 0;
        SSP.raw( 2 ).rho    = [1.5 1.5];
        Bdry.Bot.HS.rho= 0;
        SSP.raw( 2 ).alphaI = [0.5 0.5];
        Bdry.Bot.HS.alphaI = 0;
        SSP.raw( 2 ).betaI  = [ 0   0 ];
        Bdry.Bot.HS.betaI= 0 ;  
        Bdry.Top.Opt   = 'CVW .'; %S-cubic spline interpolation, C-C-linear interpolation, N-N2-linear interpolation. 
        Bdry.Bot.Opt   = 'V'; %V vacuum bellow
        Bdry.Bot.sigma = params_env{end};    
        RMax        = 0;
        Pos.s.z = sourcedepth;
        Pos.r.z = 0:deltassp:Receiverfin;      
        Beam = 'dummy';      
        write_env( envfil, model, TitleEnv, freqT, SSP, Bdry, Pos, Beam, cInt, RMax )
    end
    
    flag_env=0;
end
%% write flp

while flag_flp==1
    
    wgs84 = wgs84Ellipsoid; % reference ellipsoid
    coords_km=zeros(nnode,2);
    coords_or_km = [params_flp{5} params_flp{7}];
    for i = 1:length(lon)
        [arclen,az] = distance(coords_or(2),coords_or(1),lat(i),lon(i),wgs84); % arclen in m, az in deg
        arclen = arclen/1e3; % convert to km
        % convert azimuth angles to normal angles
        theta = 450 - az;
        if theta>360 theta = theta-360; end
        % get the coordinates in km
        coords_km(i,:) = [arclen*cosd(theta) arclen*sind(theta)]-coords_or_km;
    end
    
    fid = fopen(outfiname,'w');
    outname = sprintf("'%s'",outfiname) ;
    fprintf(fid,"'%s' 				! TITLE\n" + ...
        "'%sFM'  				! OPT\n" + ...
        "%i  	    		! NUMBER OF MODES\n" + ...
        "%i                      ! Nsx number of source coordinates in x\n" + ...
        "%.4f /  			    ! x coordinate of source (km)\n" + ...
        "%i                      ! Nsy number of source coordinates in y\n" + ...
        "%.4f /  				! y coordinate of source (km)\n" + ...
        "%i                      ! NSz\n" + ...
        "%.1f, /				! Sz( 1 : NSz ) (m)\n" + ...
        "%i                    ! NRz\n" + ...
        "%.1f %.1f /		! Rz( 1 : NRz ) (m)\n" + ...
        "%i                ! NRr\n" + ...
        "%.1f %.1f / 		! Rr( 1 : NRr ) (km)\n" + ...
        "%i                   ! Ntheta\n" + ...
        "%.1f %.1f /              ! theta( 1 : Ntheta) (degrees)\n", ...
        params_flp{1}, params_flp{2}, [params_flp{3:end-6}]);
    
    fprintf(fid,'%d					! Number of nodes\n', nnode);
    for i=1:nnode
        fprintf(fid,["%5.5f, %5.5f, '%s_%06d'\n"], coords_km(i,1),coords_km(i,2),fname_env,i);
    end
    
    m=0;
    
    fprintf(fid,'%d, 					NUMBER OF ELEMENTS (TRIANGLES)\n',length(tri(:,1)));
    for i=1:length(tri(:,1))
        m=m+1;
        if rem(i,5)~=0
            fprintf(fid,['%d, %d, %d,\n'], tri(i,1),tri(i,2),tri(i,3));
        else
            fprintf(fid,['%d, %d, %d, #%d\n'], tri(i,1),tri(i,2),tri(i,3),i);  
        end
    end
    
    fprintf( fid, "%.2f %.2f %i          /   ! ALPHA1  ALPHA2  NALPHA\n" + ... 
                "%.3f %i          /   ! STEP  NSTEPS\n" + ... 
                "%.3f          /   ! EPMULT", ...
                [params_flp{end-5:end}]);
    fclose(fid);
    
    flag_flp=0;
end

return