function p = evalri( phiS, phi, R, k, Opt, Comp )

% conversion from Fortran eval.f90
% mbp 4/2009
%
% Computes pressure field from modes
% Normalized to pressure of point source at 1 meter
% Returns the pressure matrix on a rectangular grid for a single source depth
%
% Opt = X     Cartesian   (X, z) coordinates
% Opt = R     Cylindrical (R, z) coordinates
% Opt = S     Scaled cylindrical (same as cylindrical but without the sqrt( r ) decay in range
%
% R is the range in meters

% The incoherent modal sum has not been implemented ...
% Receiver range offsets (for array tilt) have not been implemented

% For vertical component of displacement field, take finite difference in depth
% This is done using a backward difference.
% A centered difference or an FFT formula would be better

if ( Comp == 'V' )
   phidiff = diff( phi );   % should divide by diff( z ) also
   phi( 1: end - 1, : ) = phidiff;
   phi( end,        : ) = zeros( 1, size( phi, 2 ) ); % no derivative for last row
end

%phi = phi * diag( phiS, 0 );	% scale modes by phiS
phi = scalecol( double( phi ), phiS );
% form pressure field

% avoid singularity at R=0 by replacing it by a small number
Rt = R;
Rt( Rt == 0.0 ) = 1e-9;

m = length( k );
n = length( Rt );

switch Opt( 1 : 1 )
   case 'R'
      phase = spdiags( 1.0 ./ sqrt( k ), 0 , m, m ) * exp( -1i * k * Rt' ) * spdiags( realsqrt( 2 * pi ./ Rt ), 0, n, n );
   case 'X'
      phase =                                         exp( -1i * k * Rt' ) * realsqrt( 2 * pi );
   case 'S'
      phase = spdiags( 1.0 ./ sqrt( k ), 0 , m, m ) * exp( -1i * k * Rt' ) * realsqrt( 2 * pi );
end

% for horizontal component take derivative in range direction
% The following formula approximates the derivative, assuming the phase term
% (e^(-i k rr )) dominates
if ( Comp == 'H' )
   phase = diag( -1i .* sqrt( k ) ) * exp( -1i * k * Rt' ) * diag( realsqrt( 2 * pi ./ Rt ) );
end

p = phi * phase;
