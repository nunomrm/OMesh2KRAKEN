function krakelM( filename )

% Normal modes for ocean acoustics problems
% Based on the earlier Fortran version and Matlab version of Scooter
% MBP Sep. 12, 2006

clear global
tic

if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global omega SSP Bdry
global NFirstAcoustic NLastAcoustic h

% filenames
envfil = [ filename '.env' ];   % input environmental file
trcfil = [ filename '.trc' ];   % input top reflection coefficient
brcfil = [ filename '.brc' ];   % input bottom reflection coefficient

[ PlotTitle, freq, SSP, Bdry, ~, ~, cInt, ~, fid ] = read_env( envfil, 'KRAKEL' );    % read in the environmental file
fclose( fid );                  % close out the envfil

PlotTitle = ['KRAKEL(M) -' PlotTitle ];

readrc( trcfil, brcfil, Bdry.Top.Opt( 2 : 2 ), Bdry.Bot.Opt( 1 : 1 ) );   % READ Reflection Coefficients (top and bottom)

omega  = 2 * pi * freq;

% Set up vector of wavenumber samples
kMin = omega / cInt.High;
kMax = omega / cInt.Low;

NMedia          = SSP.NMedia;
h( 1 : NMedia ) = ( SSP.depth( 2 : NMedia + 1 ) - SSP.depth( 1 : NMedia ) ) ./ SSP.N( 1 : NMedia );   % vector of mesh widths
NptsAll         = sum( SSP.N( 1 : NMedia ) ) + NMedia;   % number of solution points

[ B1, B2, B3, B4, rho, Mater ] = init_elmatrix( NptsAll );    % Initialize matrices

NptsAcoustic = sum( SSP.N( NFirstAcoustic : NLastAcoustic ) ) + 1;   % size of matrix for acoustic part
[ Asparse, Bsparse ] = setAB( B1, B2, B3, B4, rho, Mater );

% note: sparse routine 'eigs' requires Bspace sym. p.d. It isn't.
%[ psi, k2 ] = eigs( Asparse, Bsparse, round( length( Asparse )/5 ), 'LR' );   %'LR' for largest real

[ psi, k2 ] = eig( full( Asparse ), full( Bsparse ) );

M = length( k2 );
%M=10

figure; imagesc( abs( psi ) );

% normalize each mode
for Mode = 1 : M
  [ z, psi( :, Mode ) ] = normiz( psi( :, Mode ), k2( Mode, Mode ), B3, B4, rho, Mater );    % normalize
end

%Min_LOC = MINLOC( Extrap( 1, 1 : M ), Extrap( 1 : M ) > omega2 / CHigh^ 2 )
%M       = Min_LOC( 1 )

k = sqrt( diag( k2 ) );

% remove any NaN's or infinities
ii = find( ~isnan( k ) );
kk = k( ii );
ii = find( ~isinf( kk ) );
k  = kk( ii );

M = length( k );

[ksort, ii ] = sort( real( k.^2 ), 1, 'descend' );
k   = k( ii );
psi = psi( :, ii );

% keep only eigenvalues with non-pos. imag. part
ii  = find( imag( k ) <= 0 );
k   = k( ii );
psi = psi( :, ii );
M   = length( k );

% Write eigenvalues to PRTFil and MODFil
fprintf( '\n\n     I          k             alpha          Phase Speed       Group Speed' )

for Mode = 1 : M
    fprintf( '\n %5i %16.10f %16.10f %16.10f', Mode, real( k( Mode ) ), imag( k( Mode ) ), omega / real( k( Mode ) ) )
end

fprintf( '\n\n' )
toc

pltitl = 'KRAKEL --';
k = k( 1 : M );
psi = psi( 1 :  length( z ), 1 : M );

save MODFIL pltitl freq k z psi Bdry NMedia rho 

field( 'MODFIL.mat' );
%movefile( 'SHDFIL.mat', shdfil );
