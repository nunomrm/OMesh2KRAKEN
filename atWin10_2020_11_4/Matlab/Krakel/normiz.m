function [ z, psiS ] = normiz( psi, X, B3, B4, rho, Mater )

% Convert the eigenvector to (U, W, TAUZX, TAUZZ) and normalize

global NMedia N
global omega
global Loc h

ik = i * sqrt( X );
NMat = length( psi );
omega2 = omega^2;

% Scale down the mode

SqNorm2 = norm( psi, 1 );
psi = psi / SqNorm2;

% Loop to compute norm of eigenvector

SqNorm = 0.0;
J    = 1;
KK   = 1;
psiS = zeros( size( psi ) );

z( 1 ) = 0.0;    % assumes DepthT = 0 for now

for Medium = 1 : NMedia
    rhoM = rho( Loc( Medium ) + 1 );
    for II = 1 : N( Medium ) + 1
        if ( strcmp( Mater( Medium, : ), 'ELASTIC ' ) ) % elastic layer
            R1 = psi( KK   );
            R2 = psi( KK + 1 );
            R3 = psi( KK + 2 );
            R4 = psi( KK + 3 );

            Contrib = h( Medium ) * ( -X * B3( J ) * R1^2 + B4( J ) * R1 * R4 - R2 * R3 );

            psiS( KK     ) = ik * psi( KK   );
            psiS( KK + 1 ) =      psi( KK + 1 );
            psiS( KK + 2 ) = ik * psi( KK + 2 );
            psiS( KK + 3 ) =      psi( KK + 3 );

            u(     II, Medium ) = ik * psi( KK   );
            v(     II, Medium ) =      psi( KK+1 );
            tauzx( II, Medium ) = ik * psi( KK+2 );
            tauzz( II, Medium ) =      psi( KK+3 );

            KK = KK+4;
        else  % acoustic layer
            Contrib = -h(Medium) * psi( KK )^2 / ( rhoM * omega2 );
            psiS( KK ) = -psi( KK );
            p( II, Medium ) = -psi( KK );
            KK = KK + 1;
            z( KK ) = z( KK - 1 ) + h( Medium );
        end

        if ( II == 1 || II == N(Medium)+1 )
            Contrib = 0.5 * Contrib;
        end
        SqNorm = SqNorm + Contrib;
        J = J+1;
    end
end

% Bottom half-space contribution

%if ( BOTOPT(1:1) == 'A' )
%    if ( strcmp( char( Mater( NMedia+1, : ) ), 'ELASTIC ' ) )
%        'Elastic halfspace normalization not implemented'
%    else
%         GAMP = sqrt( X - omega2 / CPB^2 );
%         Contrib = -psi( NMat )^2 / ( 2.0 * GAMP * rhoB * omega2 );
%     end
% end
%
%SqNorm = SqNorm + Contrib;

% Scale the mode

RN = -omega2 * SqNorm;
ScaleFac = 1.0 / sqrt( RN );
psiS( 1 : NMat ) = ScaleFac * psiS( 1 : NMat );

save rvector
