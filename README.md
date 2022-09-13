# OMesh2KRAKEN

OMesh2KRAKEN: A tool for mesh generation with OceanMesh2D for underwater acoustic modeling with KRAKEN

version 1.0 September 2022

OMesh2Kraken is under development by Nuno Monteiro (nunom@ua.pt) and Tiago Oliveira (toliveira@ua.pt). OMesh2KRAKEN is distributed under the GNU Public License.

## Dependencies
 * MATLAB (only the R2021b version was tested)
 * OceanMesh2D v4.0
 * Acoustics Toolbox (2020_11_4 version)
 * m_map package

## General description of directories
 * atWin10_2020_11_4: Contains the Acoustics Toolbox (Windows version at [1]). The original ```tests/``` folder from the Acoustics Toolbox was removed to save space;
 * data: Contains environmental data (temperature and salinity), bathymetry and shoreline datasets to run the two tests;
 * m_map: A mapping package for Matlab (dependency for OceanMesh2D) [2];
 * OceanMesh2D: OceanMesh2D toolbox [3] adapted to OMesh2Kraken. Specific adaptations were made to OceanMesh2D: the ```msh.m``` script in ```OceanMesh2d/@msh/``` was changed, and two function scripts (```writekraken3d.m```, ```interp_ss.m```) were created in ```OceanMesh2D/@msh/private```. Some folders from the original ```Examples/```, ;
 * tests: Tagus Estuary and Azores study cases shown in Sections 3.1 and 3.2 of the paper, respectivelly, to test OMesh2KRAKEN;
 * utils: Contains functions for pre/post-processing purposes (e.g., to calculate the sound speed, make bathymetry field plots), where functions ```brewermap.m``` and ```sndspd.m``` come from [4] and [5], respectively, while remaining (```find_nearest_point.m```, ```plot_mesh_bathy.m```, ```interp_ss.m``` and ```plot2D_TL_slice.m```) are custom made for OMesh2Kraken.

## How to use OMesh2KRAKEN
Here we present hands-on instructions to reproduce the Tagus estuary study case results:
1. Go to the ```tests/tagus_estuary/``` folder from the main directory and execute ```Run_OM2D_Tagus_Estuary.m```;
2. In the same directory, execute the ```Run_KRAKEN3D_Tagus_Estuary.m``` program (```KRAKEN``` and ```FIELD3D``` are run);
3. Run ```Plotting_Tagus_Estuary.m```
    
This worflow is similar for the Azores study case (```tests/azores/```), except that only ```KRAKEN``` is necessary to execute.

## Data description
 * ```data/bathymetry/EMODNET_DTM_*.nc```: bathymetry data from the EMODNET DTM database (after extraction, those two datasets were cropped for the Azores and Tagus estuary regions), that can be accessed at http://www.emodnet-bathymetry.eu;
 * ```data/bathymetry/gebco_azores_2020.nc```: bathymetry data from the GEBCO 2020 global grid dataset (https://www.gebco.net/data_and_products/gridded_bathymetry_data/gebco_2020/), used only for the Azores study case;
 * ```data/shorelines/GSHHS_h_L1*```: GSHHG shoreline data for obtained from https://www.soest.hawaii.edu/pwessel/gshhg/;
 * ```data/shorelines/CNTR_RG_01M_2020_4326*```: shoreline dataset from the GISCO database (https://ec.europa.eu/eurostat/web/gisco/overview), for the Tagus estuary study case;
 * ```data/temperature_salinity/cmems_mod_ibi_phy_anfc_0.027deg-3D_P1D-m_1647891936780.nc```: dataset of monthly mean data of 2021 of ocean's temperature and salinity 3D data in the Azores region from the Operational Mercator global ocean analysis and forecast system (https://doi.org/10.48670/moi-00016);
 * ```data/temperature_salinity/global-analysis-forecast-phy-001-024-monthly_1656015476392.nc```: 7-day span daily mean temperature and salinity data provided by the operational IBI (Iberian Biscay Irish) Ocean Analysis and Forecasting system data which was extracted for the Tagus Estuary region (https://doi.org/10.48670/moi-00027).

## References

[1] The Acoustics Toolbox. Last accessed July 8, 2022: http://oalib.hlsresearch.com/AcousticsToolbox/.

[2] Pawlowicz, R., 2020. "M_Map: A mapping package for MATLAB", version 1.4m, [Computer software], available online at www.eoas.ubc.ca/~rich/map.html.

[3] Roberts, K. J., Pringle, W. J., and Westerink, J. J., 2019. OceanMesh2D 1.0: MATLAB-based software for two-dimensional unstructured mesh generation in coastal ocean modeling, Geoscientific Model Development, 12, 1847-1868. https://doi.org/10.5194/gmd-12-1847-2019.

[4] Stephen23 (2022). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.3), GitHub. Retrieved July 8, 2022.

[5] Sergei Koptenko (2022). Sound Speed in Sea Water (https://www.mathworks.com/matlabcentral/fileexchange/4940-sound-speed-in-sea-water), MATLAB Central File Exchange. Retrieved July 8, 2022.
