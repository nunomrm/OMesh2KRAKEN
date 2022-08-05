function write_fieldsflp( flpfil, Pos )

% Write a field-parameters file

if ( ~strcmp( flpfil( end - 3 : end ), '.flp' ) )
  flpfil = [ flpfil '.flp' ]; % append extension
end

fid = fopen( flpfil, 'wt' );

rMin = min( Pos.r.range );
rMax = max( Pos.r.range );
NR   = length( Pos.r.range );

fprintf( fid, '''RP'' \t \t ! Option' );
fprintf( fid, '%6.2f %6.2f %5i \t ! rMin rMax (km) NR', rMin, rMax, NR );

fclose( fid );

