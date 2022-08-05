filename = 'tl.grid';

rd = 0 : 2 : 1000;
rr = 0 : 5 : 4000;

Nrd = length( rd );
Nrr = length( rr );

fid = fopen( filename, 'rb' );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end

TL = fread( fid, [ Nrd, Nrr ], 'float32' );    %Read complex data

fclose( fid );
whos
figure;
pcolor( rr, rd, TL );
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
shading flat;
colormap( flipud( jet ) );
colorbar
caxis( [ 40 90 ] )
