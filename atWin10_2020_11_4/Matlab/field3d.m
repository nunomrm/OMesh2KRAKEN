function field3d( filename )

% run the field3d program
%
% usage: field3d( filename )
% where filename is the field3d flp file

runfield3d = which( 'field3d.exe' );

if ( isempty( runfield3d ) )
   error( 'field3d.exe not found in your Matlab path' )
else
   eval( [ '! "' runfield3d '" ' filename ' > field3d.prt' ] )
end