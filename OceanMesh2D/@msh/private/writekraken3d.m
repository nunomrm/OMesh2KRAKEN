function outname = writekraken3d(fnames,tri,VX,B,inp_env,inp_flp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                    %
% tri: triangle (nelem,3);                                  %
% VX(:,1): longitude (nnode,1);                             %
% VX(:,2): latitude (nnode,1);                              %
% B: depth (nnode,1);                                       %
% fnames: mesh name                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nnode=length(VX(:,1));
lon=VX(:,1);
lat=VX(:,2);

flag_env = inp_env.flag;
flag_flp = inp_flp.flag;
fname_env = fnames{1};
fname_flp = fnames{2};

if flag_env==1

    disp('Writing ENV files')

    params_env = inp_env.params;
    TS_data = inp_env.TS_data;
     
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
    N=nnode;
    Nz=length(z_TS);
    Nx_d=length(lo_d);
    Ny_d=length(la_d);
    S_n=zeros(N,Nz);
    T_n=zeros(N,Nz);
    D_n=zeros(N,Nz);
    SSP_n=zeros(N,Nz);

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
    
    mesh.lon = lon;
    mesh.lat = lat;
    mesh.z = B;
    mesh.tri = tri;
    [SSP_fin,z_fin,T_fin,S_fin] = interp_ss(TS_data, mesh, deltassp);

    for idx=1:nnode
        envfil   = sprintf('%s_%06d',fname_env,idx);
        TitleEnv = envfil;
        zi=[0:deltassp:B(idx)];
        len_z=length(zi);
        SSP_idx=SSP_fin(idx,1:len_z);
        
        %Water

        if round(10.^2*zi(len_z))/10.^2~=round(10.^2*B(idx))/10.^2
            zi(end+1)=B(idx);
            SSP_idx(end+1)=[SSP_idx(len_z)];
        end

        if len_z==1
            zi(2)=deltassp;
            SSP_idx(2)=SSP_idx(1);
        end
        
        if isnan(SSP_idx(2))==1
            SSP_idx(2)=SSP_idx(1);
        end
        WalphaR=SSP_idx;     
        SSP.depth  = [0 zi(end) Receiverfin];        
        

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
        Bdry.Bot.HS.alphaR  = 0;
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
    
end
%% write flp

if flag_flp==1 

    disp('Writing FLP files')
    
    params_flp = inp_flp.params;
    coord_or = inp_flp.coord_or;
    
    wgs84 = wgs84Ellipsoid; % reference ellipsoid
    coord_km=zeros(nnode,2);
    coord_or_km = [params_flp{5} params_flp{7}];
    for i = 1:length(lon)
        [arclen,az] = distance(coord_or(2),coord_or(1),lat(i),lon(i),wgs84); % arclen in m, az in deg
        arclen = arclen/1e3; % convert to km
        % convert azimuth angles to normal angles
        theta = 450 - az;
        if theta>360 theta = theta-360; end
        % get the coordinates in km
        coord_km(i,:) = [arclen*cosd(theta) arclen*sind(theta)]-coord_or_km;
    end
    
    fid = fopen([fname_flp '.flp'],'w');

    outname = sprintf("'%s'",fname_flp) ;
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
        fprintf(fid,["%5.5f, %5.5f, '%s_%06d'\n"], coord_km(i,1),coord_km(i,2),fname_env,i);
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
end

return