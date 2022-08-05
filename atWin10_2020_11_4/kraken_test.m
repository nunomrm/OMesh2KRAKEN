clear all, close all, clc
% addpath(genpath('../windows_executables'));
addpath(genpath('windows-bin-20201102'));
addpath(genpath('input_files/'));
tic;
kraken('kraken_1500hz_50m')
elapsedTime = toc
plotshd( 'kraken_1500hz_50m.shd.mat')


%%
runkraken = which( 'kraken.exe' );
runfield  = which( 'field.exe' );

% run kraken for each of the profiles
eval( [ '! "' runkraken '" kraken_1500hz_50m' ] );