% This script runs KRAKEN and FIELD3D (KRAKEN3D) for the Tagus Estuary test case of OMesh2KRAKEN
%
% (c) 2022 Nuno Monteiro and Tiago Oliveira, University of Aveiro
%
%

clear all, close all, clc

addpath(genpath('..\..\atWin10_2020_11_4\windows-bin-20201102')); % addpath to call Windows 10 binaries of Acoustics Toolbox
fname_flp = 'tagus_estuary'; % filename of FLP file

copyfile('data_kraken/*env','.')

%%%% Run KRAKEN3D %%%%%%%%%%%

t0=tic;    % Start counting elapsed time of simulation

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
if ( isempty( runfield3d ) )
   error( 'field3d.exe not found in your Matlab path' )
else
   eval( [ '! "' runfield3d '" ' fname_flp ' > field3d.prt' ] )
end

elapsedtime=toc(t0) % Display measured elapsed time in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of run KRAKEN3D %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deletion/moving KRAKEN data files
delete *.env % ENV files are already archived at './data_kraken'
movefile('*.mod','data_kraken\');
movefile('tagus_estuary*.prt','data_kraken\');
