#SCOOTER and SPARC wavenumber integration models for tonal and broadband sources in ocean acoustic waveguides.

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

 4/92  The linear system solver (FACTOR/BACSUB) used by  SCOOTER was not initializing one of the variables to 0.  This caused sporadic NaN/Infinity problems when using IEEE arithmetic.  A similar problem occured  with the axis length specification for plots in FIELDS.

 9/92  The option of putting an acousto-elastic halfspace on the top was not working properly. (It took the density from the lower halfspace.)

10/92
FIELDS produced erroneous results when the number of wavenumber points exceeded the dimension limit of 10000.  A test has been added to flag this error and the dimension has been increased to 50000.

12/92
SPARC modified so that if the parameters are such that the time integration is fully explicit (ALPHA=BETA=0.0) then the linear solver is bypassed.  This gives some improvement in performance.

 1/93
FIELDS was writing an extra vector at zero range in cases where subsampling was required.  This caused PLOTSLICE to have problems plotting TL.

 2/93
Error in PLOTTS fixed.  The FLOW and FHIGH numbers were always used for bandpass filtering even when the option 'N' (no filtering) was selected.  If the  bandpass filter were wider than the signal bandwidth you effectively disabled filtering.  If, however,  you had a narrow bandpass filter then it was always applied to the signal and might have produced confusing results.

 3/93
Error in SOURCE fixed.  When the time-reversal option was used, SOURCE reversed the time coordinates as many times as there were sources.  As a result, the time series was reversed but the times associated with the data were not.

 6/93
Error in FIELDS fixed. Occasionally it would reference an array index that was one point beyond the last calculated value. This had no significant effect on the TL field but caused a floating point exception on machines that trap references to an undefined value.

Error in SCOOTER fixed. SCOOTER was extracting the density at the source incorrectly: it took a density at a nearby point. When the nearby point was outside the grid, an undefined density was returned.

12/93
FIELDS modified so that the specified ranges of receivers are not changed.

 3/95
The plot option in FIELDS has been removed to make SCOOTER structurally similar to KRAKEN. Similarly, the input file to fields should now be called fields.flp rather than fields.plp.

 6/95
In SPARC the option for producing a time-series on a horizontal array was not working properly. The array variable RTSTRR was being loaded instead of RTSRR.

 1/98
FIELDS automatically increases the minimum range, RMIN, to ensure that we don't try to calculate a field value right on the source. Previously RMIN was bumped to 1 m which caused problems for ultrasonic applications. Now it is bumped to RMAX / 100000. 

 7/98
FIELDS was missing the declaration of PC as a complex variable. This affected  field calculations only if the polynomial interpolation (rather than the default Pade) was used.

 7/98
A test has been added to SPARC so that it will abort if the user specifies a larger mesh than SPARC has been dimensioned for.

3/3/01
The Pade interpolation option (actually the default) has occasionally produced NaN's. Following the modification in Numerical Recipes, a small perturbation has been introduced to avoid the 'rare 0/0 condition'. The problem appears to have been eliminated.

5/01
I needed the top and bottom reflection coefficient options in SCOOTER. This had not been exercised in many years (if ever) and several problems showed up. Open statements for the files were missing and the section of code that does dynamic allocation for the top reflection coefficient was using the variables for the bottom reflection coefficient, so it failed if they were already in use. Finally the DOS script needed modifying.

6/01
Looks like the option for using a pre-calculated reflection coefficient in SCOOTER had fallen into disrepair (or was never originally tested). Fixes were made to make the read process consistent with the current format of the IRC file used for that.

12/03
Some clean up in the Hankel transform routine (HTS.f90). Also, the option to select dB vs. linear scale in FIELDS.f90 has been eliminated since FIELDS.f90 no longer plots anything directly. Checks have been added during the read-in stage for that first option line so that a user using an old input file will get a clear error message. Finally, FIELDS has been implemented in MATLAB (see at/Matlab/fields2.m). Since FIELDS appears to be a reserved word in MATLAB the script is currently fields.m. SCOOTER and SPARC used to fold each row of the Green's function matrix across records. This was done to 1) keep the record length below the 4094 maximum in old f77 compilers, and 2) to keep the memory storage down for that same matrix by buffering it to the output file as each record was filled. Since this is a nuisance; since the f95 does not appear to have a maximum record length; since real memory is cheap and virtual memory seems like it ought to be fairly efficient here, I've eliminated that. Finally, the format of the Green's function file (GRNFIL) has been made identical to that of a shade file (SHDFIL).

5/05
Variable 'IniFlag' had not been declared as type LOGICAL in fields.f90. The defaul INTEGER was fine, except for the G95 compiler, which detected use of an integer as a logical ... fixed. Also, the subroutine FTS.f90 (inside HTS.f90) had fallen into disrepair. This is used for the rarely-used line-source option in FIELDS.f90. fixed.

6/05
The main tridiagonal matrix in SCOOTER has been converted to double precision. Roundoff errors were causing some problems when the number of elements was made large. (The problem showed up using 5000 points over 100 m in a Pekeris waveguide with a 50 Hz source. Thus we were using over 1000 points per wavelength, which you'd probably not normally do. However, you should be able to do this without breaking things.)

5/06
Changes to FIELDS to bring it in line with the Matlab version.

1/08
The Pade interpolation routine is used in FIELDS to interpolate the Green's function kernel. I had found that it seemed to act up in some cases and made Polynomial interpolation the default. However, my logic for making that the default was off, so the Pade method was being used (caused a problem for at/tests/halfspace when the wavenumber sampling was reduced. The logic has been fixed.

11/13
Added an option to zero out the stabilizing attenuation. This is useful for noise calculations.

February 2017
Testing with an aeroacoustic scenario (involving low density) revealed that there was a missing density factor in handling tabulated reflection coefficients. This occured in both KRAKEN and SCOOTER for this option and has been fixed.

March 2017
The formula for picking the index of the source density was wrong, causing an error when the density was depth dependent and not equal to 1. This showed up in an atmospheric modeling case where the source was in air and with a wild depth dependence. The fix that has been implemented might still have a minor offset if the source is not in the first medium.

November 2017
fieldsco.m has been modified to allow an irregularly spaced vector of receiver ranges. This is convenient for noise modeling. This required a change to the input structure, which also now matches the style used for all the other models in the Acoustics Toolbox.

An option 'H' has also been added to use the exact Bessel transform, which is important in the very near field.

January 2018
The broadband implementation in SCOOTER had not been completed. Many changes were made to make that operational.

December 2018
A user had a problem running fields.f90 with a broadband SCOOTER run. In fact, fields.f90 had not been modified to handle that. Code was added to have fields.f90 terminate with an error when this is attempted.

A bug was found in this latest version of fields.f90 in that it was not allocating storage for the source x and y coordinates. Scooter doesn't use those variables except to write them to the output file to have a consistent format with the 3D models that do use that information. This had no real impact; however, the Intel Fortran compiler detected it when all the debug checks were enabled.

February 2018
In converting the codes to broadband, I introduced a bug in fields.f90 where the frequency variable, Freq, was not initialized. This broke the option for applying a source beam pattern. Note that fields.f90 is not really supported. fieldsco.m (the Matlab equivalent) is the one that normally should be used.

May 2019
The routine sourceMod.f90 is used by sparc.f90 to read a source timeseries file. A line copied a matrix of these timeseries into a single vector for subsequent filtering. The vector was longer than the corresponding entry in the matrix. With runtime diagnostics turned on in the Intel compiler this triggered an error about accessing data outside the array bounds. Those values beyond the array bound were never used. However, this was fixed.

The Matlab routine fieldsco.m has been modified to allow negative values in the range vector. This doesn't generally make physical sense, but there are occasions where we like to show such symmetric plots.

June 2019
For a broadband run, SCOOTER does automatic scaling of the finite-element grid based on the frequency. An error was fixed where the scaling was done using the lowest frequency as the reference, even though the original grid was set up at the central frequency given in the second line of the env file.

May 2020
A user specified that there was one receiver depth, but then provided two values in the following line. This caused confusion in the Matlab version of Scooter in terms of the vector allocated for the receiver depths. The routine readvector.m was modified to ignore extraneous extra depths.

July 2020
The fieldsco.m routine is the Matlab routine usually used to transform the scooter Green's function to the pressure field. Normally the far-field approximation to the Hankel function was used. In that case, there is a singularity at the origin. If the user selected a zero (or small) range it was shifted to avoid that singularity. However, that shift was also being applied when a user selected the exact Hankel function. So for near-field problems, with that option, an error was introduced. Fixed ....

August 2020
(See also March 2017.) A quick fix had been done in March 2017 to allow a source in the atmosphere over the ocean. The issue here was to correct for the density the source. The fix did not allow the possibility of sources below the first medium, if it had a different density. This caused errors in a particular test case. Code has been added to more carefully identify the density at the source depth. However, we note that there is an ambiguity about which density to use if sources are placed on interfaces between media with different densities.

October 2020
Originally the Acoustics Toolbox considered fields that depended on source depth (Sz), receiver depth (Rz), and receiver range (Rr), i.e. a 3D array. Changes were required for BELLHOP3D which allows a ( Sx, Sy, Sz ). Also a broadband option was added bringing in a timestep. This has all been unified in the Fortran codes.

The Matlab codes have not been kept current with all those capabilities and sparcM was not in synch with the scooterM format. sparcM has been updated and is now consistent with scooterM. However, further work will be needed for broadband scooterM runs.
