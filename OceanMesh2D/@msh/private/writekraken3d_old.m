function outname = writekraken3d(outfiname,tri,VX,B,title,coords_or,calc,params_flp,fname_env,TS_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                    %
% tri: triangle (nelem,3);                                 %
% VX(:,1): longitude (nnode,1);                             %
% VX(:,2): latitude (nnode,1);                              %
% B: depth (nnode,1);                                       %
% outfiname: mesh name                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('write FLP mesh and ENV files');

nnode=length(VX(:,1));
lon=VX(:,1);
lat=VX(:,2);

fname_env = 'tagus_estuary';

% T-S data extraction

lo_d=TS_data.lon;
la_d=TS_data.lat;
T=TS_data.T;
S=TS_data.S;
z_TS=TS_data.z;
x1=min(lon); x2=max(lon); y1=min(lat); y2=max(lat);
x_f=lo_d>=x1&lo_d<=x2;
y_f=la_d>=y1&la_d<=y2;
lo_d=lo_d(x_f);
la_d=la_d(y_f);


% Calculation of SSP in the unstructured mesh from T-S data

N=nnode;
Nz=length(z_TS);
Nx_d=length(lo_d);
Ny_d=length(la_d);
S_n=zeros(N,Nz);
T_n=zeros(N,Nz);
D_n=zeros(N,Nz);

% figure, set(gcf,'Units','normalized','outerposition',[.01 .05 .8 .75])
for i_z=1:Nz
    S_z=squeeze(S(x_f,y_f,i_z,1));
    T_z=squeeze(T(x_f,y_f,i_z,1));
    for i_m=1:N
        diff=zeros(Nx_d,Ny_d);
        for i = 1:Nx_d
            for j = 1:Ny_d
                diff(i,j)=abs(lo_d(i)-lon(i_m))+abs(la_d(j)-lat(i_m));
            end
        end
        [idx_x, idx_y]=find(diff==min(diff(:)));
        c=0;
        while 1
            if isnan(S_z(idx_x,idx_y)) 
                if c==numel(S_z)-1
                    S_n(:,i_z)=ones(N,1).*NaN;
                    T_n(:,i_z)=ones(N,1).*NaN;
                    D_n(:,i_z)=ones(N,1).*NaN;
                    break
                else
                    diff(idx_x,idx_y)=1e10;
                    [idx_x, idx_y]=find(diff==min(diff(:)));
                    if length(idx_x)>1  idx_x=idx_x(1);  end
                    if length(idx_y)>1  idx_y=idx_y(1);  end
                    S_n(i_m,i_z)=S_z(idx_x,idx_y);
                    T_n(i_m,i_z)=T_z(idx_x,idx_y);
                    D_n(i_m,i_z)=z_TS(i_z);
                    c=c+1;
                end
            else
                S_n(i_m,i_z)=S_z(idx_x,idx_y);
                T_n(i_m,i_z)=T_z(idx_x,idx_y);
                D_n(i_m,i_z)=z_TS(i_z);
                break
            end
        end
    end
    Si = squeeze(S(:,:,i_z,1));
%     subplot(121)
%     hold off
%     pcolor(lo,la,Si')
%     colormap jet
%     hold on
%         caxis([30 max(S(:))]); %c=colorbar;
% 
%     x = [x1, x2, x2, x1, x1];
%     y = [y1, y1, y2, y2, y1];
%     plot(x, y, 'k-', 'LineWidth', 3);
%     plot(m.p(:,1),m.p(:,2),'k.')
%     xlim([x1-0.05 x2+0.05])
%     ylim([y1-0.05 y2+0.05])
%     subplot(122)
%     hold off
%     trisurf(tri,lon,lat,S_n(:,i_z)')
%     hold on, shading interp, grid off
%     trimesh(tri,lon,lat,ones(N,1).*1e6, 'facecolor', 'none', 'edgecolor', 'k')
%     caxis([30 max(S(:))]); c=colorbar;
%     c.Label.String = 'Salinity';
%     colormap(jet)
%     colorbar
%     xlim([x1-0.05 x2+0.05])
%     ylim([y1-0.05 y2+0.05])
%     view(0,90)
%     sgtitle(['level=' num2str(i_z)])
%     pause(3)
end

% remove NaN layers
c=0;
for i=1:Nz
    if sum(isnan(T_n(:,i-c)))==N
        T_n(:,i-c)=[];
        S_n(:,i-c)=[];
        D_n(:,i-c)=[];
        c=c+1;
    end
end

% sounspeed profiles
SSP_n=sndspd(S_n,T_n,D_n);


D=B; % topo-bathymetry
disp('Writing .env files')
ir=0;
cInt.Low = 1400;
cInt.High = 15000;
model    = 'KRAKENC';
Receiverfin = round(max(D)+0.1*max(D)); 
deltassp=1;
Receivernumber=500;
sourcedepth=15;
for idx=1:N
  envfil   = sprintf('%s_%06d',fname_env,idx);
  TitleEnv = envfil;
  freqT    = 1500;

  SSP.NMedia = 2;
  SSP.N      = [ 10000 10000 ];
  SSP.sigma  = [ 0   0 ]; %Hrms

  
  %Water
  if D(idx)>=deltassp
    zi=0:deltassp:D(idx); 
    if round(10.^2*zi(length(zi)))/10.^2==round(10.^2*D(idx))/10.^2
        ftest(idx)=1;   
    else
        ftest(idx)=0;
        zi(length(zi)+1)=D(idx);
    end
  elseif D(idx)<deltassp
    if D(idx)<2e-2
        zi=[1e-2 2e-2];
    else
        zi=[1e-2 D(idx)];
    end
  end
  SSP.depth  = [ 1e-2 zi(end) Receiverfin ];
    %   zi(length(z)+1)=D(idx);
  
    %WalphaR=ones(size(zi,2),1)*1500; %interp1(SSP_modelKraken.z(1,:),reshape(ssp(yy,i_m,:),1,length(SSP_modelKraken.z)),zi);      
    WalphaR=interp1(D_n(idx,:),SSP_n(idx,:),zi,'spline');
    f_z=find(zi>max(D_n(:)));
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
    SSP.raw( 2 ).alphaR = [ 1700 1700 ]; % speed
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
    Bdry.Bot.sigma = 1700;      
    RMax        = 0;
    Pos.s.z = sourcedepth;
    Pos.r.z = linspace( 0, Receiverfin, Receivernumber );      
    Beam = 'dummy';      
    write_env( envfil, model, TitleEnv, freqT, SSP, Bdry, Pos, Beam, cInt, RMax )
end


%% write flp

wgs84 = wgs84Ellipsoid; % reference ellipsoid
coords_km=zeros(nnode,2);
coords_or_km = [params_flp(3) params_flp(5)];
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
disp( title )  ;
fprintf(fid,"'%s' 				! TITLE\n", title);
fprintf(fid, "'%sFM'  				! OPT\n" + ...
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
    calc,params_flp(1:end-6));

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
            params_flp(end-5:end));
fclose(fid);
return