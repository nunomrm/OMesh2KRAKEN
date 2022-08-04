function outname = writeflp(outfiname,EToV,VX,B,title,fname_env,coords_km,calc,props)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:                                                    %
% EToV: triangle (nelem,3);                                 %
% VX(:,1): longitude (nnode,1);                             %
% VX(:,2): latitude (nnode,1);                              %
% B: depth (nnode,1);                                       %
% outfiname: mesh name                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('write FLP mesh');

node(:,1)=(1:length(VX(:,1)));
node(:,2)=VX(:,1);
node(:,3)=VX(:,2);
node(:,4)=B;
title
fid = fopen(outfiname,'w');
outname = sprintf("'%s'",outfiname) ;
disp( title )  ;
fprintf(fid,"'%s' 				! TITLE\n", title);
fprintf(fid, "'%s'  				! OPT\n" + ...
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
    calc,props(1:end-6));

fprintf(fid,'%d					! Number of nodes\n', length(node(:,1)));
for i=1:length(node(:,1))
%     fprintf(fid,['%5.5f, %5.5f, %s_%06d.env\n'], VX(i,1),VX(i,2),fname_env,i);
    fprintf(fid,["%5.5f, %5.5f, '%s_%06d'\n"], coords_km(i,1),coords_km(i,2),fname_env,i);
end

m=0;

fprintf(fid,'%d, 					NUMBER OF ELEMENTS (TRIANGLES)\n',length(EToV(:,1)));
for i=1:length(EToV(:,1))
    m=m+1;
    if rem(i,5)~=0
        fprintf(fid,['%d, %d, %d,\n'], EToV(i,1),EToV(i,2),EToV(i,3));
    else
        fprintf(fid,['%d, %d, %d, #%d\n'], EToV(i,1),EToV(i,2),EToV(i,3),i);  
    end
end

fprintf( fid, "%.2f %.2f %i          /   ! ALPHA1  ALPHA2  NALPHA\n" + ... 
            "%.3f %i          /   ! STEP  NSTEPS\n" + ... 
            "%.3f          /   ! EPMULT", ...
            props(end-5:end));
fclose(fid);
return