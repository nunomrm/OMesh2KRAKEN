MODULE KrakenMod

  USE MathConstants
  SAVE

  INTEGER,            PARAMETER :: ENVFile = 5, PRTFile = 6, MODFile = 20, EVMFile = 22, MaxMedium = 500, NSets = 5
  LOGICAL                       :: CountModes
  INTEGER                       :: FirstAcoustic, LastAcoustic, NV( NSets ), ISet, M, &
                                   LRecordLength, IRecProfile = 1, ModeCount, Mode, IProf, ifreq
  REAL    (KIND=8)              :: ET( NSets ), hV( NSets ), cMin, cLow, cHigh, freq, omega, omega2, RMax
  REAL    (KIND=8), ALLOCATABLE :: EVMat( :, : ), Extrap( :, : ), VG( : )
  COMPLEX (KIND=8), ALLOCATABLE :: k( : )
  CHARACTER (LEN= 8)            :: TopOpt, BotOpt
  CHARACTER (LEN=80)            :: Title

  ! finite-difference grid
  INTEGER                       :: Loc( MaxMedium ), NG( MaxMedium ), N( MaxMedium )
  REAL      (KIND=8)            :: h(   MaxMedium )

  ! storage for finite-difference equations
  REAL (KIND=8),    ALLOCATABLE :: B1( : ), B1C( : ), B2( : ), B3( : ), B4( : ), rho( : )

END MODULE KrakenMod
