function outname = writekraken3d(tri,VX,B,inp_env,inp_flp)
%
% (c) 2022 Nuno Monteiro, University of Aveiro
%
% Writes required input ENV and FLP files for KRAKEN and 
% FIELD3D, from an OceanMesh2D unstructured mesh and given
% temperature and salinity data for sound speed calculation
% (the latter not required, and a constant value of sound
% speed can be attributed to all ENV files)
%
%%% INPUTS
% tri - triangular mesh elements from the mesh, with N_e 
% elements: (N_e,3)
% VX - mesh vertices (or nodes) coordinates (nnode,2)
% B - bathymetry field (nnode, 1)
% inp_env - struct with inputs for ENV file creation
% inp_flp - struct with inputs for FLP file creation
%
%%% OUTPUTS:
% outname - output filename
%

% number of nodes, longitude and latitude
nnode=length(VX(:,1));
lon=VX(:,1);
lat=VX(:,2);

% extract flags for ENV/FLP file writing
flag_env = inp_env.flag; % 0 - write ENV files, 1 - don't write ENV files
flag_flp = inp_flp.flag; % 0 - write FLP file, 1 - don't write FLP file

if flag_env == 1

    disp('Writing ENV files')
    
    % extract struct fields from inp_env 
    fname_env = inp_env.fname;
    params_env = inp_env.params;
    deltassp=params_env.deltassp;
    z_fin = 0:deltassp:max(B);
    cInt.Low = params_env.C_lim(1);
    cInt.High = params_env.C_lim(2);
    model = 'KRAKEN';
    Receiverfin = params_env.Rdepth; 
    sourcedepth = params_env.Sdepth;
    freqT = params_env.freq;
    SSP.NMedia = params_env.NMedia;
    SSP.N =  params_env.N_layers;
    SSP.sigma = params_env.sigma_int;
    
    % Conditions for sound speed calculation or not (depending on presence of temperature/salinity input data)
    if isfield(inp_env,'TS_data')==1 % activation of sound speed calculation
        TS_data = inp_env.TS_data;
        % Calculation of SSP in the unstructured mesh from T-S data
        mesh.lon = lon;
        mesh.lat = lat;
        mesh.z = B;
        mesh.tri = tri;
        [SSP_fin, T, S] = calc_ss(TS_data, mesh, z_fin);
        
    else    % in case a sound speed constant value for all points (inp_env.SS_constant) or nothing is given
        if isfield(inp_env,'SS_constant')==1
            SSP_fin = ones(nnode,length(z_fin)).*inp_env.SS_constant;
        else
            SSP_fin = ones(nnode,length(z_fin)).*1500;
        end
    end
    
    % Water sound speed profiles (SSPs): adjustments
    for idx = 1:nnode
        envfil = sprintf('%s_%06d',fname_env,idx);
        TitleEnv = envfil;
        
        zi = 0:deltassp:B(idx);
        len_z = length(zi);
        SSP_idx = SSP_fin(idx,1:len_z);        
        if zi(len_z) ~= B(idx) && len_z>1 % condition to insert the bathymetry at each node as the final depth and attribute 
                                        % the sound speed value of the deepest point
            zi(end+1) = B(idx);
            SSP_idx(end+1) = SSP_idx(len_z);
        elseif len_z == 1 && zi < deltassp  % condition to assure two depth levels in the water columns of ENV files
                                            % (if there is one depth, KRAKEN returns error)
            if zi(1) > 1e-2
                zi = [0 B(idx)];
                SSP_idx = [SSP_idx(1) SSP_idx(1)];
            else
                zi = [0 1e-2];
                SSP_idx = [SSP_idx(1) SSP_idx(1)];
            end
        end
        
        if str2num(sprintf('%.2f',zi(end)))==str2num(sprintf('%.2f',zi(end-1)))
            zi(end) = [];
            SSP_idx(end) = [];
        end
        
        % create adequate struct files to save ENV files with the help of write_env (function from the Acoustic's Toolbox)
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
        SSP.raw( 2 ).alphaR = params_env.C_bot;
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
        Bdry.Bot.sigma = params_env.sndspd_bot;    
        RMax        = 0;
        Pos.s.z = sourcedepth;
        Pos.r.z = 0:deltassp:Receiverfin;      
        Beam = 'dummy';      
        write_env( envfil, model, TitleEnv, freqT, SSP, Bdry, Pos, Beam, cInt, RMax )
    end
    
end

if flag_flp==1 

    disp('Writing the FLP file')
    
    % load data from inp_flp struct
    fname_flp = inp_flp.fname;
    params_flp = inp_flp.params;
    coord_or = inp_flp.coord_or;
    
    coord_km=zeros(nnode,2);
    coord_or_km = params_flp.coord_or;
    
    % conversion of mesh node coordinates from lon/lat (degrees) to km
    for i = 1:length(lon)
        [arclen,az] = distance(coord_or(1,2),coord_or(1,1),lat(i),lon(i),wgs84Ellipsoid); % arclen in m, az in deg
        arclen = arclen/1e3; % convert to km
        % convert azimuth angles to normal angles
        theta = 450 - az;
        if theta>360 theta = theta-360; end
        % get the coordinates in km
        coord_km(i,:) = [arclen*cosd(theta) arclen*sind(theta)]-coord_or_km;
    end
    
    % write to FLP file
    fid = fopen([fname_flp '.flp'],'w');
    outname = sprintf("'%s'",fname_flp);
    fprintf(fid,"'%s'    ! TITLE\n" + ...
        "'%sFM' ! OPT\n" + ...
        "%i ! NUMBER OF MODES\n" + ...
        "%i ! Nsx number of source coordinates in x\n" + ...
        "%.2f %.2f / ! x coordinate of source (km)\n" + ...
        "%i ! Nsy number of source coordinates in y\n" + ...
        "%.2f %.2f / ! y coordinate of source (km)\n" + ...
        "%i    ! NSz\n" + ...
        "%.1f %.1f / ! Sz( 1 : NSz ) (m)\n" + ...
        "%i ! NRz\n" + ...
        "%.1f %.1f / ! Rz( 1 : NRz ) (m)\n" + ...
        "%i ! NRr\n" + ...
        "%.2f %.2f / ! Rr( 1 : NRr ) (km)\n" + ...
        "%i ! Ntheta\n" + ...
        "%.1f %.1f / ! theta(1 : Ntheta) (degrees)\n", ...
        params_flp.title, params_flp.calc, params_flp.Nm, ...
        params_flp.Ns(1), params_flp.Sxy_lim(1,1), ...
        params_flp.Sxy_lim(1,2), params_flp.Ns(2), ...
        params_flp.Sxy_lim(2,1), params_flp.Sxy_lim(2,2), ... %
        params_flp.NSz, params_flp.Sz_lim(1), ...
        params_flp.Sz_lim(2), params_flp.NRz, ...
        params_flp.Rz_lim(1), params_flp.Rz_lim(2), ...
        params_flp.NRr, params_flp.Rr_lim(1), ...
        params_flp.Rr_lim(2), params_flp.NRtheta, ...
        params_flp.Rtheta_lim(1), params_flp.Rtheta_lim(2));
    
    fprintf(fid,'%d    ! NUMBER OF NODES\n', nnode);
    for i=1:nnode
        fprintf(fid,"%5.5f, %5.5f, '%s_%06d'\n", coord_km(i,1),coord_km(i,2),fname_env,i);
    end
    m=0;
    fprintf(fid,'%d !    NUMBER OF ELEMENTS (TRIANGLES)\n',length(tri(:,1)));
    for i=1:length(tri(:,1))
        m=m+1;
        if rem(i,5)~=0
            fprintf(fid,'%d %d %d\n', tri(i,1),tri(i,2),tri(i,3));
        else
            fprintf(fid,'%d %d %d #%d\n', tri(i,1),tri(i,2),tri(i,3),i);  
        end
    end
    
    fprintf( fid, "%.2f %.2f %i          /    ! ALPHA1  ALPHA2  NALPHA\n" + ... 
                "%.3f %i          /    !    STEP  NSTEPS\n" + ... 
                "%.3f          /    ! EPMULT", ...
                params_flp.GBtheta_lim(1), params_flp.GBtheta_lim(2), ...
                params_flp.NGBtheta, params_flp.GBstep, ...
                params_flp.GBsteps, params_flp.epsilon_mult);

    fclose(fid);
end

return
