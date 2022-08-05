#KRAKEN Normal modes for ocean acoustics

   Copyright (C) 2009 Michael B. Porter

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


UPDATES:

 1/30/86  New versions of KRAKEN and KRAKENC correct a bug in
          eigenvector computations which probably did not affect
          you unless you saw the error message "failure to
          converge in SINVIT" 

 2/18/86  Bug in KRAKENC corrected which affected eigenvector
          accuracy for strongly attenuated modes.  When
          significant, this error would have been flagged as
          "SINVIT FAILED TO CONVERGE". The field computation has
          been made a common subroutine called EVAL for use by
          PLOTTLR, PLOTTLD, FIELD, ... to assure consistency and
          simplify any changes to the field computation. The
          modified version stores the field matrix in a fashion
          which allows for up to 50000 field points essentially
          independent of whether the matrix is 50000x1 or
          1x50000, etc. With EVAL a common routine, the X,R or S
          parameter which indicates cartesian, cylindrical or
          scaled cylindrical field computations has been
          incorporated in all the field computing routines. 

 3/12/86  TILT and DISPLACE merged into the FIELD program so that
          one routine is used for vertical, tilted or displaced
          arrays. Also, a blank line for last record in '.RED',
          '.SOD' and '.RER' files is no longer needed. 

 4/8/86   CLOW parameter added for compatibility with SCOOTER FFP
          program. At some later date KRAKEN will be modified to
          allow the mode calculation to be restricted to the
          [CLOW,CHIGH] phase speed interval, but for now CLOW is
          a dummy parameter. 

 4/23/86  KRAKEN and KRAKENC modified to allow one or more
          elastic media above the water column.  Thus an elastic
          ice model may be used. Additionally 'KRAKEN- ' or
          'KRAKENC- ' is now inserted before the user supplied
          title.  With several programs now producing shade files
          it is desireable to have this indicator on plotted
          output. 

 5/22/86  FIELD and PLOTTLR modified to treat range dependent
          environments by adiabatic mode theory. 

 6/1/86   Surface scatter options modified to accomodate Twersky
          ice scatter as an amplitude only effect.  Hard surface
          backing option removed since it seemed more a hazard
          (of accidental use) than a benefit.  Some unit numbers
          changed: most input is now done from unit 5 rather than
          1 which makes it simpler to run the program
          interactively since unit=5 is the default input unit;
          the unit numbers for '.VEC' and '.WNO' have been
          changed so that with adiabatic mode runs a series of
          such files can be allocated to sequential unit numbers. 

 6/15/86  FIELD3D module added which generalizes adiabatic mode
          theory to 3-dimensional problems. 

 11/1/86  Environmental file closed promptly after data is read for
          user convenience. START module added to provide a start-up
          field for PE programs. 

 11/30/86 Spline option fixed.  This had apparently not been used
          by anybody and was not properly updated when the input
          file format was changed.  Also, added check to
          gracefully abort the program when the problem is too
          big for dimensioned storage. 

 1/8/87   Slight changes in format statements on input to avoid
          the problem of neglected decimal points causing
          disasters. 

 2/13/87  KRAKENC modified to gracefully handle problems where
          the user requests more modes than KRAKENC is
          dimensioned to handle. Both KRAKEN and KRAKENC modified
          to narrow a user specified phase speed interval which
          is unnecessarily large. (Eliminates the problem of a
          user specifying a CLOW which is so low as to cause a
          failure to converge in the root finder.) Dimensions
          increased to allow for up to 1000 modes. 

 3/20/87  All modules modified to use a single mode file instead
          of separate files for eigenvalues and eigenvectors. 
          Thus '.WNO' and '.VEC' are now '.MOD'.  This simplifies
          the command files, especially for 3D runs which employ
          numerous mode files. Adiabatic versions of PLOTTLR and
          FIELD have been made the only versions.  Thus PLOTTLRAD
          replaces PLOTTLR.  This is a step in the direction of
          reducing the proliferation of new modules with minor
          variations over existing modules. SCOOTER FFP program
          integrated into the package. 

 6/1/87   Options for different kinds of boundary conditions
          added. Allows free, rigid and homogeneous half-spaces
          at top and bottom. Environmental file format changed
          slightly to accomodate same. 

 8/10/87  Bug in section which scales to avoid under/over flow
          corrected. Caused an overflow error to occur in the
          special case where a mode became exactly 0.0 at one of
          the finite difference points. 

 2/1/88   Errors in line-source represention corrected. The line
          source is of virtually no relevance in underwater
          acoustics so this probably doesn't affect you.  Coupled
          mode solution implemented but density variations and
          half-space contributions are ignored. 

 6/29/88  Inefficiencies in mode reading removed. Source,
          receiver depth files eliminated.  ( Use of the /
          command to terminate records combined with an
          interpolation routine made the '.SOD' and '.RED' files
          an unnecessary incumberance. ) KRAKEN input structure
          modified so that all lines are now read in using list
          directed I/O. 

 8/1/88   Some problems in tabulated bottom refl. coeff.
          corrected. Also corrections to FIELDS which made it
          necessary to use more k-space points than it should
          have needed. Interfacial scatter finally implemented
          correctly. 

 10/5/88  Modified to allow up to 301 SSP pts.  Changes to
          internal structures to do this efficiently.   

 6/10/89  Addition of source/receiver depth line.  This allows
          the user to specify the depths at which modes are
          written to the mode file.  Changes to various plotting
          routines to allow the user to specify the axis lengths. 
          Change to SSP input structure: in the old version the
          user had to provide a count of the number of SSP points
          in each layer.  In the new version you simply provide
          the depth at the bottom of each layer.  The program
          keeps reading SSP points till it reaches that depth.

 9/19/89  Several structural changes to the codes have been made.
          A subroutine READIN replaces the previous INCLUDE for
          reading the environmental file and a routine BCIMP has
          been written which combines TOPIMP and BOTIMP routines
          for the boundary conditions.  BCIMP has been written so
          that there are no differences between the options
          available for top and bottom boundaries. (Actually,
          that's not entirely true since the internal reflection
          coefficient option is still only implemented for the
          lower boundary.)  The interfacial roughness calculation
          has been separated from the mode normalization code to
          improve readability. 

 10/30/89 The option 'M' for attenuation in dB/m has been
          implemented. This is a one-line formula in the routine
          CRCI which, despite appearing in earlier documentation,
          seems to have gone unimplemented till now.

 11/20/89 Changes in FIELD3D calculations.  In the original
          version the receiver was located at the origin and the
          source moved out in range.  The present version is
          consistent with the other field calculation routines in
          that the source is fixed at the origin and the
          receivers move out in range.  In addition, the option
          for multiple source depths is now implemented.

 12/11/89 KRAKEN/KRAKENC modified to produce a MODFIL even when
          there are no propagating modes.  This eliminates the
          need to introduce DUMMY nodes in a FIELD3D run.

 2/2/90   Changes to mode file format. The new
          format is much more compact for problems with lots of
          modes. In addition, it includes the density and
          half-space properties so that the field can be
          calculated in a half-space.  Also FIELD now allows a
          title to be specified.

          Error in FIELD removed. In range-dependent calculations
          the number of modes, M, may sometimes be reduced in
          propagating out in range.  The error in FIELD was to
          not reset M for subsequent sources, possibly causing
          erroneous results.  This error would only have affected
          range-dependent problems with multiple sources.

          Changes to KRAKEN and KRAKENC to reduce storage.

          Small change in adiabatic calculations (routines EVALAD
          and EVAL3D.  The denominator in the adiabatic
          expression involves sqrt(k(r)r).  In the previous
          version (as in certain other modal codes the mean value
          of k(r) was used instead.  In practice, you'll be
          hard pressed to spot the difference in TL calculations.

 3/30/90  Errors in PROFIL routine corrected.  The LOC variable
          was not being saved between calls and there was also an
          illegal reference to LOC(0).  On a VAX the former
          problem is a non-issue because VMS Fortran
          automatically saves variables.  The latter is only a
          problem when LOC(0) accesses a non-zero value or a
          protected memory location.

          Also, if you ever used the Thorp attenuation option you
          know that it wasn't working right.  A missing factor
          was causing attenuation values about 10,000 times too
          large!

 6/20/90  Bug in FIELD3D corrected which caused errors when the
          step size between points where TL was calculated was
          large enough that there was not at least one TL point
          in each triangle. An analogous problem was also
          corrected in FIELD for adiabatic calculations.  PDQ
          option added to FIELD3D for quicker but less accurate
          3D calculations.  Finally, FIELD3D has been modified to
          allow radials to pass directly over nodes of the
          triangulation.  Similarly, the source can now lie
          anywhere in the plane (on a node or even outside the
          triangles).

 12/10/90 Extensions to coupled mode routine.  Density and halfspace
          contributions are now included.  This has necessitated
          some modifications to the MODFIL format.

          A program BOUNCE has been created to compute a reflection
          coefficient for a multilayered domain.  This reflection
          coefficient can be passed to KRAKENC to provide a boundary
          condition or plotted using the new routine PLOTRTH.

          NEW IMPROVED KRAKENC! Up to 3 times slower than the old
          version!  I've added an option which causes KRAKENC to
          perform a very careful root search.  In certain elastic
          problems the old root finder was found to miss modes.
          Unfortunately making the root finder more robust entails
          slowing the code down.  As a compromise, I've introduced
          the slow/robust root-finder as an option.  For more
          information see TOPOPT as described in KRAKEN.HLP.

  1/26/91 The XL and XR arrays were erroneously dimensioned one
          element too small.  This caused an error in the special
          case when the number of modes was exactly equal
          to MAXM, the maximum number of modes KRAKEN is dimensioned
          to calculate.

          Modifications to the adiabatic mode calculation routine
          EVALAD.  With a very coarse distribution of receivers in
          range the results differed from the convergent answer
          obtained with say 500 range points.  Basically, there was
          some sloppy handling of the phase integral in range
          which assumed you had a dense sampling of points in
          range.  (Practically the inaccuracy introduced was
          probably less than that due to the adiabatic
          approximation itself.)


  5/15/91 Error corrected in READIN.  The vectors CP and CS were
          dimensioned to a length of one element.  When the analytic
          profile option is used this is not sufficient.

  5/20/92 Common blocks reshuffled to avoid alignment errors on
          various UNIX systems during compilation.  Double opens
          implemented to get record lengths from files on UNIX
          systems.  A missing EXTERNAL statement has been added to
          the SLATECBESSEL package.  On SGI systems this was flagged
          at compilation time and subsequently led to erratic
          results.

          An error has been corrected in KRAKEN which
          caused minor discrepancies with the KRAKENC result.
          KRAKEN was taking the single precision part of the SSP
          rather than the full double precision available.  Since
          the models write the mode files in single precision, this
          error only shows up on the print-out where the full precision
          is shown.  Also, it did not affect results on VAX machines.

          An error in the 'Q' type attenuation has been fixed.  The
          error caused any attenuation to be ignored when using this
          option.

 9/25/92  The option of putting an acousto-elastic halfspace
          on the top was not working properly. (It took the
          density from the lower halfspace.)

 11/9/92  Error corrected in coupled-mode calculations (subroutine
          EVALCM).  The program multiplied the profile ranges
          by a factor of 1000 to convert km to m every time it
          was called. The error was that the scaling was repeated
          for sources after the first.  This had the effect of
          stretching the environment so much that it was effectively
          range-independent for subsequent sources.

 1/20/93  Form feeds removed.  TYPE and ACCEPT statements removed.
          Routine COVAR written to generate a covariance matrix
          for matched-field processing with colored or white
          noise. BART and CAPON implemented to compute
          ambiguity surfaces.  A few other minor changes were
          made to accomodate the more sensitive AIX (IBM) and
          AbSoft (NeXT) compilers.

 6/ 5/93  KRAKEN and KRAKENC have been modified to discard
          any eigenvalue returned by the root finder if
          an error occurred during the search. This fixes
          a problem that occurred in some elastic cases
          where extra modes were returned. Also, the e(i pi/4)
          factor has been added to the pressure calculations.
          Since this phase does not affect TL or standard
          beamformers, it had never been included.

 8/20/93  KRAKENC has been modified to zero out the imaginary
          part of an eigenvalue if it's positive.  A positive
          imaginary part may occur when the true root is just
          below the real axis. For such roots the root-finder
          can give a very accurate approximation to the root
          that just happens to be on the positive side of
          the axis.

          A check has been added in COVAR for zero eigenvalues.
          The colored noise model will always fail in such
          cases since without attenuation the noise sheet
          leads to infinite intensities.

12/26/93  Unix version was not creating empty mode files
          correctly for the case when no modes were present.

 4/ 9/94  Adiabatic calculations were showing excessive
          round-off error when many (>1000) receiver
          ranges were used.  The offending variable has
          been converted to double precision in EVALAD.F and
          EVALD3D.F

          KRAKENC was not setting the value for the
          last calculated mode when it aborted the
          run because the mesh was too coarse.

          All programs have been set up to open Fortran files
          by name rather than unit number. This gets around
          the problem that different vendors use different
          default names for files causing the scripts to
          need to be modified for each machine.

 6/23/94  Bug fix in KRAKEN: a variable was not properly initialized
          in the special case where a mode cutoff on the 3rd
          or subsequent meshes. To fix it, FUNCT has been modified
          to return 0 and zsecx.f has been modified to quit immediately
          when f(x) = 0 (identically).

 3/13/95  Various routines for calculating the acoustic pressure
          were missing a factor of i (eval.f, eval3d.f, evalcm.f, ...
          This has no effect on TL but leads to a Hilbert transform
          in time-series calculations.

 3/15/96  The logical record length used in opening the mode files
          was too small in the case of many layers and few receivers.

 4/19/96  getmod.f was testing the proximity of the receiver depth
          to the tabulation grid using points above and below the
          receiver. This is not valid when there is only one tabulation
          point ... The code has been modified to bypass linear
          interpolation for this case.

 5/27/96  An option has been added to calculate the finite difference
          grid (NMESH) automatically.

10/29/96  For adiabatic field calculations, erroneous results are
          obtained if the first profile does not start at zero range.
          A test has been added to field.f to enforce this.

10/14/97  FIELD3D failed without explanation in a case where there were
          more modes than allocated storage (MAXM). A test has been added
          to flag this error in all routines.

11/30/97  MODBIN/MODASC have been modified so that the ascii file format
          is easier to read into matlab programs. LAPACK implementation
          of RAN() included to aid portability. Bug in BOUNCE ...
          it was accessing (but not using) LOC( 0 ) in calculating the
          reflection coefficient for a stack of elastic layers with no
          acoustic layer.

 1/98 PLOTRTH has been modified to include new display options.

 8/98 KRAKEN and KRAKENC now allow up to 1000 profiles to be stacked one after
      another in a single ENVFIL. Thus, a coupled-mode or adiabatic run can be
      done with a single call to KRAKEN or KRAKENC. In many cases, the run-time
      was dominated by just the cost of loading KRAKEN so this can be a big time
      saver ...

 8/98 The inverse iteration routines (sinvitd and sinvitz for KRAKEN and KRAKENC
      respectively, have been cleaned up and tuned a bit for about a 25% speed
      improvement.

 3/%& Group velocity is now included in the output print file.
      Modifications for Y2k compliance.

 6/00 The option for having KRAKEN and KRAKENC do multiple runs in one envfil was not working (as a result of my having added a 'close' for the file in the wrong place (fixed). Variables also needed to be deallocated so that an error was not generated in trying to allocate an already exiting variable. (This problem crept in when features of f90 dynamic allocation were introduced.)
 
10/00 Evidently some errors crept into FIELD3D/EVAL3D when dynamic array allocation was added. The range-vector had not been allocated. Subroutine 'out' was carrying the vector 'iset' through to a calling program (not using it itself) and had never allocated it as a vector. There was an erroneous attemp to allocate vector 'sumk' for subsequent radials, when it had already been allocated.

10/00 The documentation files have been updated to reflect the minor changes in the input structure required for dynamic allocation.

2/01 Bottom and top reflection coeffiecient files were not being closed. This causes a problem when multiple profiles were stacked in a single KRAKEN envfil and the reflection coefficients needed to be re-read for subsequent profiles.

9/01 Same problem as on 6/00 but with modfil's not being properly closed. As I understand the Fortran spec, all files should be automatically and properly closed on program termination. However, since the compiler doesn't seem to be doing this properly I'm adding close statements for all files.

11/01 A bug in Bounce was introduced in converting to f95. The bug involved the calculation of the mesh spacing H(). Fixed. An additional feature has been added allowing BOUNCE to directly output a reflection coefficient file in the '.brc' format that BELLHOP uses.

May 2002
Attenuation units can be input in terms of a 'loss parameter'. Incomplete dependencies in Makefiles have been corrected.

June 2002
An atan(y/x) operation was converted to atan2(y, x) to allow for x=0. This happens when KRAKEN or SCOOTER, calculateS a reflection coefficient for vertical incidence.

November 2002
Changes to KRAKEN so that mode files are numbered modfil0001, modfil0002, etc. as opposed to modfil1, modfil2, ... This addresses a problem that the (presumed standard) I0 format option in Lahey Fortran does not exist in Portland Group Fortran. 

December 2002
The field program allows you to specify receiver displacements (for tilted arrays). Typically, the array is straight and one wants to just put in a value of zero and have it repeated for all phones. The FIELD code was inadvertantly accessing an unitialized variable in that case. I've never seen a problem result from this, but it's obviously not the way to leave things and it was being trapped when I had all the debug traps set by the compiler. It's now fixed.

February 2003
A problem was brought to my attention in using the option for a tabulated reflection coefficient in KRAKENC. It turned out that the code was using the sound speed at the top of the ocean, rather than the value at the water/sediment interface. The phase of the reflection coefficient was also being included with an incorrect sign (relative to the convention in BOUNCE). The first of these problems was also present in SCOOTER (but not the second). A set of test problems (at/tests/TabRefCoef) has been added to exercise this option. There are 6 cases there; however, BELLHOP is not currently part of that battery, nor are the top reflection coefficients tested ... KRAKENC and SCOOTER are tested for both BRC and IRC types of reflectioncoefficients.

March 2003
In an earlier process of converting the code to f95 I had used a cmplx( ,,kind=4) statement where kind=8 (double precision) should have been used. (This is in function CRCI.) The effect of this was that the sound speed profiles were provided in single precision. Given that in the real world we know the sound speed that precisely, the error is irrelevant for practical problems. However, it would effect extremely precise benchmarking exercises.

May 2004
Compilation under the latest version of Lahey Fortran revealed problems related to not deallocating dynamic variables in subroutines that were called more than once. DEALLOC statements have been added to KRAKEN and KRAKENC in the necessary places.

May 2006
The dynamic allocation done in FIELD3D to take advantage of fortran95 had some small problems. In the process of looking at this, a large number of changes were done. The shade files now allow multiple bearings which affects all the codes in the Acoustics Toolbox. An ElementMod (Element Module) has been introduced to take further advantage of the f95 dynamic allocation. EVALGB.f90 (the Gaussian beam option for horizontal refraction) has had key sections converted to double precision since a test case (at/test/3DAtlantic) showed sensitivity to the precision.

November 2006
The default mesh spacing of 10 points/wavelength was found too coarse for calculating the reflection coefficient using BOUNCE. It was increase to 20 points/wavelength. Unfortunately this effects all the models in the toolbox that use that same readin routine. Those models will now run slower than necessary.

December 2007
Gulf test case using the adiabatic option generates a segmentation fault with the Intel compiler. g95 and gfortran are fine. No errors are detected by the Intel compiler when all checks are turned on. Error is eliminated by using fewer receiver depths. Speculation is that the automatic array allocation for the variable PhiL is failing (since this is one difference compared to the coupled mode version).

April 2008
Fixed logic that was supposed to detect the last SSP point by an approximate match to the value of Depth. The test was too stringent, and didn't always find a match because of roundoff.

Learned that the Intel compiler generates a segmentation fault when the stack overflows. A compiler switch was added to the Makefile to replace the stack by a heap for the Intel compiler. This answers issue of December 2007 above ...

Added logic to gracefully exit when an end of file is detected in the ENVFIL. This happens, for instance, when a user puts extra blank lines at the end. KRAKEN thinks the user has concatenated runs and tries to process the blank lines as an input stream.

March 2009
The "subtab" option for profiles in field.f90 had become broken ... This option allows you to specify just the first and last ranges of the profiles in an adiabatic or coupled-mode calculation. Then field would subtabulate intermediate ranges. The call to the subtab routine was just in the wrong place ...

August 2009
Previously, KRAKEN and KRAKENC would create a different mode file for each profile in the environmental file. The new version writes the modes into a single, combined mode file. All the routines that read mode files, including FIELD have been modified accordingly.

The mode files use a fixed record length determined by the first profile in the environmental profile. Therefore an error will be generated if the required record length increases in subsequent profiles. The record length is normally controlled by the number of distinct source/receiver depths.

The Matlab versions of FIELD have also been converted to be compatible with this new format. The coupled mode implementation in Matlab has been modified to store the modes in single precision (they are written to the mode file in single precision, so there is not a lot of advantage in computing later in double precision).

The convergence test in the root finders was tighter in KRAKENC than KRAKEN (root finders zseccx and zsecx, respectively). These have been unified using the tighter test, but with a looser tolerance for the test. Also, there was a reliable process in zseccx to avoid NaNs in the secant method when the function returned the same value at adjacent iterates. This has been carried over to zsecx. This should be a more robust process.

The inverse iteration routines (invited and sinvitz) were not able to get the specified growth rate for certain interfacial modes. The test has been relaxed a bit.

The test in KRAKEN and KRAKENC to determine when a mode was outside the user-specified spectral integral was incorrect for negative eigenvalues. Fixed.

Constants, such as 1.0 have been typed in as 1.0D0 everywhere to make sure they are treated as double precision.

April 2011
Modified FIELD3D to allow multiple source locations in the x-y direction  (for noise calculations with sources located everywhere).

May 2011
BOUNCE modified to handle imaginary angles more gracefully when writing a reflection coefficient file. They are mapped to zero degrees.

November 2011
The MergeVectors routine in at/misc was supposed to remove duplicates in the source/receiver depths where the modes are tabulated. It was not consistently doing so because it tested for exact equality. This has been corrected.

November 2012
The routine (getmodes.f90) that returns the mode shapes in halfspaces used a formula
that was only valid with KRAKENC modes. There had actually been a
comment in the code noting that limitation ... That routine has been
modified to detect whether the modes were calculated by KRAKEN or
KRAKENC and use the appropriate formula.

June 2015
There are both Fortran and Matlab implementations of coupled modes in
KRAKEN. The Matlab versions had a formula for calculating the
integration weights that assumed a single index for the interface was
returned. However, for internal interfaces a special weight is needed
for points on either side of the interface. This caused a fatal error
when the Matlab code was run iff there were sub-bottom layers. This
affected both evalcm.m and evalcmLoadAll.m which are two different
implementations of coupled modes.

September 2015
There is a rarely used and little tested option in FIELD3D to use Gaussian beams x modes to treat
horizontal refraction in the adiabatic mode approximation. The file involved is evaluategb.f90.

A shortcut
had been taken in the code in assuming that a single range step along
the central ray of the beam could only take you into an adjacent
element. Code has been added to have it skip through multiple elements
as necessary.

The code had also been ignoring the imaginary part of the modal wavenumbers (the attenuation factor). (The
complex rays are projected onto the real plane by taking the real part
of the wavenumbers; however, this should not have been applied to the
phase delay.)

February 2017
Testing with an aeroacoustic scenario (involving low density) revealed
that there was a missing density factor in handling tabulated
reflection coefficients. This occured in both KRAKEN and SCOOTER for
this option and has been fixed.

November 2017
A new option had been added in KRAKEN and KRAKENC some months back to
do broadband calculations. This is invoked with an option letter and
requires a frequency vector to be added to the input file. This change
was fairly involved, requiring modifications to the file formats for
both modfiles and shdfiles.

Two bugs showed up. First, the record length did not consider the
length of the frequency vector. When there were many frequencies, the
record lenght was too short to accomodate that vector.

Second, the automatic grid selection took a scalar multiple based on
the frequency. However, it was scaling things for each frequency based
on how much larger it was relative to the lowest frequency. However,
the initial grid itself was selected based on the highest
frequency. Therefore, the grids were ridiculously fine and the runtime
was excessive.

January 2018
The format of the mode files has been modified so that top and bottom
halfspace information is written for each frequency and for each
profile. The previous version assumed that did not change. The
soundspeed generally changes with frequency because it is a complex
number incorporating the frequency-dependent attenuation. Allowing the
sound speed in the halfspace also to change with profile makes KRAKEN
more general. These changes affect both the writing of the mode files
as well as Fortran and Matlab routines that read those mode files.

The KRAKEN Field routines have also been pulled out into a separate
directory to organize things a bit more cleanly.

September 2020
A user found inconsistent results between the Fortran and Matlab
version of field. This was because they had applied array tilt and the
Matlab version does not have that implemented. (The Fortran version does.) A test was added so
that that is flagged as an error.

October 2020
All the references to the Twersky ice scatter model have been
removed from both the code and the documentation. This capability had
fallen out of repair decades ago and (in my opinion)
was of minimal value.

At some point, perhaps two years ago, I broke the Matlab version of
the coupled mode code. I had added code that interpolated the ranges
of new profiles based on a first and last range. However, that was
generally not appropriate. The Fortran version remained as the
default. This problem has been fixed and both codes produces the same
results. This corrects resutls from at/tests/Gulf and at/tests/step
when using the Matlab coupled mode code.
