OMesh2KRAKEN: A tool for mesh generation with OceanMesh2D for underwater acoustic 
modeling with KRAKEN

#######################################################################################
#######################################################################################

OMesh2KRAKEN is distributed under the GNU Public License.

v.0.0 August 2022 * OMesh2Kraken is under development by Nuno Monteiro (nunom@ua.pt) 
and Tiago Oliveira (toliveira@ua.pt) at the University of Aveiro.

#######################################################################################
#######################################################################################

Dependencies:
 * MATLAB (only the R2021b version was tested)
 * OceanMesh2D v4.0
 * Acoustic's Toolbox (2020_11_4 version)
 * m_map package

#######################################################################################

General description of directories:
 * atWin10_2020_11_4: Contains the Acoustics Toolbox (Windows version at [1]). The 
original "tests" folder from the Acoustics Toolbox was removed to save space;
 * data: Contains environmental data (temperature and salinity), bathymetry and 
shoreline datasets to run the two tests;
 * m_map: A mapping package for Matlab (dependency for OceanMesh2D) [2];
 * OceanMesh2D: OceanMesh2D toolbox [3] adapted to OMesh2Kraken. Changes were 
made to msh.m in OceanMesh2d/@msh/, and a script writekraken3d.m to save input files
for KRAKEN3D was created in OceanMesh2D/@msh/private;
 * tests: Tagus Estuary and Azores study cases shown in Sections 3.1 and 3.2 of the
paper, respectivelly, to test OMesh2KRAKEN;
 * utils: Contains functions for pre/post-processing purposes (e.g., to calculate the 
sound speed, make bathymetry field plots), where functions "brewermap.m" and "sndspd.m"
come from [4] and [5], respectively, while remaining ("find_nearest_point.m, plot_mesh_
bathy.m" and "plot2D_TL_slice.m") are custom made for OMesh2Kraken.

#######################################################################################

Test step-by-step instructions (Tagus estuary study case):
 * Tagus estuary study case
    1) Go to the tests/tagus_estuary folder from the main directory and execute 
  "Run_OM2D_Tagus_Estuary.m";
    2) In the same directory, execute Run_KRAKEN3D_Tagus_Estuary.m (KRAKEN and FIELD3D
are run);
    3) Run Plotting_Tagus_Estuary.m 
This worflow is similar for the Azores study case (except that in the second step only
KRAKEN is runned).

#######################################################################################

Data description:
 * "data/bathymetry/EMODNET_DTM_*.nc": bathymetry data that was extracted from the 
EMODNET DTM database (after extraction, those two datasets were cropped for the Azores 
and Tagus estuary regions), that can be accessed at http://www.emodnet-bathymetry.eu
 * "data/shorelines/GSHHS_h_L1*": GSHHG shoreline data for obtained from 
https://www.soest.hawaii.edu/pwessel/gshhg/
 * "data/shorelines/CNTR_RG_01M_2020_4326*": shoreline dataset from the GISCO database
(https://ec.europa.eu/eurostat/web/gisco/overview), for the Tagus estuary study case
 * "data/temperature_salinity/cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m_1647891936780
.nc": dataset of monthly mean data of 2021 of ocean's temperature and salinity 3D data 
in the Azores region from the Operational Mercator global ocean analysis and forecast 
system (https://doi.org/10.48670/moi-00016)
 * "data/temperature_salinity/global-analysis-forecast-phy-001-024-monthly_16560154763
92.nc": daily 3D temperature and salinity data from the from the operational IBI 
(Iberian Biscay Irish) Ocean Analysis and Forecasting system; data was extracted for 
the Tagus Estuary region with a 7-day span (and daily frequency) during September 2020 
and comes from https://doi.org/10.48670/moi-00027

#######################################################################################

References:

[1] The Acoustics Toolbox. Last accessed July 8, 2022: http://oalib.hlsresearch.com/AcousticsToolbox/.

[2] Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

[3] Roberts, K. J., Pringle, W. J., and Westerink, J. J., 2019. OceanMesh2D 1.0: MATLAB-based software for two-dimensional unstructured mesh generation in coastal ocean modeling, Geoscientific Model Development, 12, 1847-1868. https://doi.org/10.5194/gmd-12-1847-2019.

[4] Stephen23 (2022). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub. Retrieved July 8, 2022.

[5] Sergei Koptenko (2022). Sound Speed in Sea Water (https://www.mathworks.com/matlabcentral/fileexchange/4940-sound-speed-in-sea-water), MATLAB Central File Exchange. Retrieved July 8, 2022.
