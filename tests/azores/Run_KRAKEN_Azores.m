% This script runs Acoustics Toolbox's KRAKEN model for the Azores test case of OMesh2KRAKEN
%
% (c) Nuno Monteiro and Tiago Oliveira, University of Aveiro
% September 2022
%

clear all, close all, clc

addpath(genpath('..\..\atWin10_2020_11_4\windows-bin-20201102')); % addpath to call Windows 10 binaries of Acoustics Toolbox

copyfile('data_kraken/*env','.','f')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Run KRAKEN %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

elapsedtime=toc(t0) % measured elapsed time to run kraken in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% end of run KRAKEN %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deletion/moving KRAKEN data files
delete *.env % ENV files are already archived at './data_kraken'
delete *.mod
movefile('azores*.prt','data_kraken\');
