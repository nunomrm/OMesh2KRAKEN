function [SSP_fin,z_fin,T_fin,S_fin] = interp_ss(TS_data, mesh_data, deltassp)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % T-S data extraction
    lo_d=TS_data.lon;
    la_d=TS_data.lat;
    T=TS_data.T;
    S=TS_data.S;
    z_TS=TS_data.z;

    % mesh data extraction
    lon=mesh_data.lon;
    lat=mesh_data.lat;
    tri=mesh_data.tri;
    B=mesh_data.z;
    N=length(B);

    x1=min(lon); x2=max(lon); y1=min(lat); y2=max(lat);
    x_f=lo_d>=x1-0.1&lo_d<=x2+0.1;
    y_f=la_d>=y1-0.1&la_d<=y2+0.1;
    lo_d=lo_d(x_f);
    la_d=la_d(y_f);
    [LO_D,LA_D]=meshgrid(lo_d,la_d);
    LO_D=LO_D';
    LA_D=LA_D';
    
    % Calculation of SSP in the unstructured mesh from T-S data
    Nz=length(z_TS);
    if max(z_TS)>max(B)
        f=find(z_TS>max(B));
        Nz=f(1);
        z_TS=z_TS(1:Nz);
    end
    Nx_d=length(lo_d);
    Ny_d=length(la_d);
    S_n=zeros(N,Nz);
    T_n=zeros(N,Nz);
    D_n=zeros(N,Nz);
    SSP_n=zeros(N,Nz);
    for i_z=1:Nz
        S_z=squeeze(S(x_f,y_f,i_z,1));
        T_z=squeeze(T(x_f,y_f,i_z,1));
        S_zz=zeros(Nx_d*Ny_d,1);
        T_zz=zeros(Nx_d*Ny_d,1);
        lo_zz=zeros(Nx_d*Ny_d,1);
        la_zz=zeros(Nx_d*Ny_d,1);
        c=1;
        for i=1:Nx_d
            for j=1:Ny_d
                S_zz(c)=S_z(i,j);
                T_zz(c)=T_z(i,j);
                lo_zz(c)=LO_D(i,j);
                la_zz(c)=LA_D(i,j);
                c=c+1;
            end
        end

        la_zz=double(la_zz);
        lo_zz=double(lo_zz);
        nan_flags=isnan(S_zz);
        if sum(nan_flags)<Nx_d*Ny_d
            la_zz(nan_flags)=[];
            lo_zz(nan_flags)=[];
            S_zz(nan_flags)=[];
            T_zz(nan_flags)=[];
        else
            if i_z>1
                S_n(:,i_z)=S_n(:,i_z-1);
                T_n(:,i_z)=T_n(:,i_z-1);
            elseif i_z==1
                while 1
                    S_zz=squeeze(S(x_f,y_f,i_z,1));
                    sum_nan_flags=sum(isnan(S_zz));

                    if sum_nan_flags<Nx_d*Ny_d
                        break;
                    elseif sum_nan_flags==Nx_d*Ny_d & i_z==Nz
                        error('No data in the input T-S grid!');
                        break;
                    else
                        i_z=i_z+1;
                    end
                end
            end
        end
        
        F1=scatteredInterpolant(lo_zz,la_zz,S_zz,'linear','nearest');
        F2=scatteredInterpolant(lo_zz,la_zz,T_zz,'linear','nearest');
        warning(''); % clear last warning msg
        S_int = F1(lon,lat);
        [warnMsg1, warnId1] = lastwarn;
        T_int = F2(lon,lat);
        [warnMsg2, warnId2] = lastwarn;
        if ~isempty(warnMsg1) | ~isempty(warnMsg2) % to prevent deficient interpolations detected by rare warning messages: "The underlying triangulation is empty"
            S_n(:,i_z) = S_n(:,i_z-1);
            T_n(:,i_z) = T_n(:,i_z-1);
        else
            S_n(:,i_z) = F1(lon,lat) ;
            T_n(:,i_z) = F2(lon,lat) ;
        end
        for i_m=1:N
            if B(i_m)<z_TS(i_z)
                T_n(i_m,i_z)=NaN;
                S_n(i_m,i_z)=NaN;
            end
        end
    end
    
    z_fin=[[0:deltassp:deltassp:max(B)]];
    SSP_fin = zeros(N,length(z_fin)).*NaN;
    T_fin = zeros(N,length(z_fin)).*NaN;
    S_fin = zeros(N,length(z_fin)).*NaN;

    for idx=1:N
        if isnan(T_n(idx,2))==1
            T_n(idx,2)=T_n(idx,1);
        end
        if isnan(S_n(idx,2))==1
            S_n(idx,2)=S_n(idx,1);
        end
        zi=0:deltassp:deltassp:B(idx);
        if length(zi)==1
            zi(2)=deltassp;
        end

        
        ff=find(z_TS>max(zi));
        if sum(ff)>0
            j = ff-1;
        else
            if j>1
                j = length(zi)
            else
                j = 2;
            end
            break
        end
        z_inp=z_TS(1:j)';
        T_inp=T_n(idx,1:j);
        S_inp=S_n(idx,1:j);
%         if length(z_inp)==1
%             z_inp(2)=z_TS(2);
%             T_inp(2)=T_inp(1);
%             S_inp(2)=S_inp(1);
%         end

        T_idx=interp1(z_inp,T_inp,zi','linear','extrap');
        S_idx=interp1(z_inp,S_inp,zi','linear','extrap');

        SSP_idx=sndspd(S_idx,T_idx,zi');
        T_fin(idx,1:length(zi))=T_idx;
        S_fin(idx,1:length(zi))=S_idx;
        SSP_fin(idx,1:length(zi))=SSP_idx;
        
    end
  
end