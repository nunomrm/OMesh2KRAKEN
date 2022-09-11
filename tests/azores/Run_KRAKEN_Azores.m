%% initialize and copy env files to current directory

clear all, close all, clc

addpath(genpath('..\..\atWin10_2020_11_4\windows-bin-20201102')); % to call windows model executables

copyfile('data_kraken/*env','.','f')

%% run kraken

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

elapsedtime=toc(t0) % measured elapsed time to run kraken in seconds

%% delete or move kraken input/output files
delete *.env % env files are already archived at './data_kraken'
delete *.mod
movefile('azores_*.prt','data_kraken\');