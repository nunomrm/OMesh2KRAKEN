function fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, option, Rminkm, Rmaxkm, Nrr )
%
% calculates the field using the Green's function produced by SCOOTER
%
% usage: fieldsco( filename )
% You must include the file extension, if it exists
% Reads control info from a file fields.flp
% mbp, Dec. 2003

% differences between this one and the Fortran version:
% The Fortran code uses a recursion to interpolate G(k) in the complex
% k-plane, using Pade or Polynomial fits.
% That code should map directly into Matlab with minimal effort; however,
% as a recursion it will likely run very slowly in Matlab.
%
% The compromise used here does interpolation of G(k) as if k were real
% numbers, using Matlab's interp1 routine to do piecewise linear
% interpolation.
% The main drawback there is that there is no interpolation with interp1
% with respect to the imaginary component of k. This means that we need to
% keep the sampling for the interpolated wavenumbers very close to the
% original sampling used by SCOOTER in generating G(k).
%
% Bottom line: if the user oversamples G(k) when running SCOOTER, then the
% Matlab version of FIELDSCO will also use more points in the FFT than it
% needs to.

% fft sampling is complicated ...
% keep in mind that the usual FFT has f in [0, fmax) and t in [0, tmax)
% note the open interval; r( fmax ) = r( 0 ) and r( tmax ) = r( 0 )
% following lines are user input on desired range sampling

% phase speed limits for optional tapering
% this is for user play (at your own risk)
cmin = 1e-10;
cmax = 1e30;

% if user is just supplying a file name then read all the variables
if ( nargin == 1 )
    % read fields.flp data to select range limits
    fid = fopen( 'fields.flp', 'r' );
    option = fgetl( fid );   % x or r coordinate; positive/negative/both sides of spectrum
    option = option( 2:3 ); % only two letters of option line are used
    Rminkm = fscanf( fid, '%f', 1 );
    Rmaxkm = fscanf( fid, '%f', 1 );
    Nrr    = fscanf( fid, '%i', 1 );
    fclose( fid );

    % read in the Green's function
    [ PlotTitle, ~, freq, atten, Pos, Gtemp ] = read_shd( filename );
    k = Pos.r.range;
end

Nsd = length( Pos.s.depth );    % # of source depths
Nrd = length( Pos.r.depth );
Nk  = length( k );              % # of wavenumbers

fprintf( '\n\n--- FieldSCO --- \n\nRminkm = %d, Rmaxkm = %d, Nrr = %i \n', Rminkm, Rmaxkm, Nrr )
if ( Rminkm < 0 )
    %warning( 'Rmin <0, G(-r) defined as G(r) if r<0' )
    error( 'Rmin must be nonnegative' )
end
Rmax    = 1000 * Rmaxkm;
Rmin    = 1000 * Rminkm;
deltar  = ( Rmax - Rmin ) / ( Nrr - 1 );
rr      = Rmin : deltar : Rmax;
RmaxFFT = Rmax + deltar;

kleft  = 2.0 * pi * freq / cmax; % left  limit for tapering
kright = 2.0 * pi * freq / cmin; % right limit for tapering

% fiddle with transform to satisfy user's desired range sampling
[ kInterp, Nt, deltakInterp, rr, Nrr, deltar, IRatiodeltar ] = set_transform( k, Nk, Rmin, deltar, Nrr, RmaxFFT );

% keep only the ranges in the user specified interval and subsample if deltar was bumped
NrrLast      = min( round( ( Rmax - Rmin ) / deltar )+1, Nt );
   
for isd = 1: Nsd
   G = squeeze( Gtemp( isd, :, : ) );
   if size( G, 2 ) == 1 % if G is a vector, expand it into a matrix with one row
      G = reshape( Gtemp, 1, length( Gtemp ) );
   end
  
   G = taper( G, k, Nk, kleft, kright );
   % interpolate Green's function on to this new grid
   GInterp = interp1( k, G.', kInterp, 'linear', 0 ).';  % uses value of 0 for kInterp outside the grid
   
   if ( option( 1: 1 ) == 'X' )
      ptemp = FTS( kInterp( 1 ), deltakInterp, atten, rr, GInterp, option );
   else
      ptemp = HTS( kInterp( 1 ), deltakInterp, atten, rr, GInterp, option );
   end
   p( isd, :, : )  = ptemp( :, 1 : IRatiodeltar : NrrLast );    % subsample
   fprintf( 'Transform completed for source depth %f \n', Pos.s.depth( isd ) );
end

Pos.r.range  = rr(       1 : IRatiodeltar : NrrLast );
atten = 0;

PlotType = 'rectilin  ';
if ( length( filename ) >=6 && strcmp( filename(1:6), 'GRNFIL' ) )
   save SHDFIL PlotTitle PlotType freq atten Pos p
else
   shdfilname = filename(1:end-4)
   save( shdfilname, 'PlotTitle', 'PlotType', 'freq', 'atten', 'Pos', 'p' )
end

end % end of fieldsco

%%******************************************************

function [ kInterp, Nt, deltakInterp, rr, Nrr, deltar, IRatiodeltar ] = set_transform( k, Nk, Rmin, deltar, Nrr, RmaxFFT )

%%******************************************************
% set up wavenumber sampling
%%******************************************************

% deltak is what SCOOTER used; deltakInterp is for interpolation
% If deltar is too big, take submultiple

kMax = k( 1 ) + 2.0 * pi / deltar;   % kMax required to get users specified deltar

% If that kMax is less than what was used in SCOOTER, force the user to
% sample the range more finely (later we'll subsample in range to give what
% was requested)

IRatiodeltar = 1;   % this must always be returned by set_transform
if ( kMax < k( end ) )
   IRatiodeltar = ceil( ( k( end ) - k( 1 ) ) / ( kMax - k( 1 ) ) );
   deltar = deltar / IRatiodeltar;
   Nrr    = IRatiodeltar * ( Nrr - 1 ) + Nrr;
   kMax   = k( 1 ) + 2.0 * pi / deltar;
   fprintf( 'Number of ranges, Nrr, increased to %i so that wavenumber limit exceeds kMax used by SCOOTER \n', Nrr );
end

% If necessary, increase Nt to ensure that deltak is fine enough to take
% you to the largest receiver range

Nt2 = round( RmaxFFT * ( kMax - k( 1 ) ) / ( 2 * pi ) );  % DeltaK limit
if ( Nt2 > Nrr )
   fprintf( 'Transform size, Nt, increased from NR = %i to %i. \n', Nrr, Nt2 );
   fprintf( 'Thus we are zero filling the wavenumber spectrum to effectively interpolate on to a finer grid \n' )
end
Nt = max( Nrr, Nt2 );

% bump Nt if necessary to make sure deltakInterp is not coarser than deltak grid
deltak       = ( k( end ) - k( 1 ) ) / Nk;
deltakInterp = 2 * pi / ( Nt * deltar );

if ( deltakInterp > deltak )
   IRatio = round( deltakInterp / deltak );
   Nt     = IRatio * Nt;
   fprintf( 'Transform size, Nt, bumped to %i to ensure deltak sampling is fine enough \n', Nt );
end

% find a good transform size that is at least as large a NT

for ii = 1:100  % search 100 integers for factors
    number_of_factors( ii ) = length( factor( Nt+ii-1 ) );
end

[ junk, ii ] = max( number_of_factors );    % take the one that had the most factors
Nt = Nt + ii - 1;
%Nt = 2^ceil( log2( Nt ) );  % make Nt a power of 2 for faster FFT
fprintf( 'Transform size, Nt, bumped to %i for fast FFT \n', Nt );

deltakInterp = 2.0 * pi / ( Nt * deltar );
kInterp = k( 1 ) : deltakInterp : k(1) + ( Nt - 1 ) * deltakInterp;

kMax    = kInterp( end ) + deltakInterp;  % because the fft goes from 0 to kmax but G(kmax) is not computed
Nt      = length( kInterp );
deltar  = 2 * pi / ( kMax - kInterp( 1 ) );
Nrr     = Nt;
rr      = Rmin: deltar : Rmin + ( Nrr - 1 ) * deltar;
end % of set_transform

%%******************************************************

function G = taper( G, k, Nk, kleft, kright )

% windowing to smooth out any discontinuities at the end of the spectrum

if ( kleft > k( 1 )  )
    Nwinleft = 2 * round( ( kleft - k( 1 ) ) / ( k(end) - k( 1 ) ) * Nk ) + 1;    % odd number covering about 10% of the spectrum
    winleft  = hanning( Nwinleft )';
    Nwinleft = ( Nwinleft - 1 ) / 2;
    winleft = winleft( 1: Nwinleft ); % taking just the left half of the symmetric Hanning window
else
    Nwinleft = 0;
    winleft  = [];
end

if ( kright < k( end )  )
    Nwinright = 2 * round( ( k(end) - kright ) / ( k(end) - k( 1 ) ) * Nk ) + 1;    % odd number covering about 10% of the spectrum
    winright  = hanning( Nwinright )';
    Nwinright = ( Nwinright - 1 ) / 2;
    winright  = winright( end - Nwinright + 1: end ); % taking just the right half of the symmetric Hanning window
else
    Nwinright = 0;
    winright  = [];
end

window  = [ winleft ones( 1, Nk - Nwinleft - Nwinright ) winright ];
G = scalecol( G, window );
end % of taper
