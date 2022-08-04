%% initialize _soupy env files to current directory

clear all, close all, clc

addpath(genpath('..\..\atWin10_2020_11_4\windows-bin-20201102')); % to call windows model executables

copyfile('io_files_kraken/*env','.')

%% run kraken and field

t0=tic;

runkraken = which( 'kraken.exe' );

DirInfo = dir( '*.env' );

for ifile = 1 : length( DirInfo )
    [ ~, filename, ~ ] = fileparts( DirInfo( ifile ).name );
    if ifile==1
        disp( filename )
    elseif rem(ifile,10)==0
        disp( filename )
    end
    eval( [ '! "' runkraken '" ' filename ] );
end

disp('Running FIELD3D ...')
runfield3d = which( 'field3d.exe' );
filename = 'tagus_estuary';
if ( isempty( runfield3d ) )
   error( 'field3d.exe not found in your Matlab path' )
else
   eval( [ '! "' runfield3d '" ' filename ' > field3d.prt' ] )
end

elapsedtime=toc(t0) % measured elapsed time to run kraken+field3d in seconds

%% deletion/moving input/output files of KRAKEN
delete *.env % env files are already archived at './io_files_kraken'
movefile('*.mod','io_files_kraken\');
movefile('tagus_estuary_*.prt','io_files_kraken\');