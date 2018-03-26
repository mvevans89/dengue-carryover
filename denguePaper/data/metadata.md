# Metadata for Data for Evans et al. 'Carry-over effects'

Note that all `.RData` files were saved via `saveRDS`, so must be read in as `readRDS`, **not** `load`.

Data begins in raw folder and is transformed via `01_climate_processing.Rmd` or `02_mosquito_processing.Rmd` into clean data. Some values are calculated as noted in description.

## emergence

### raw

`augustEmergence.csv`|`octoberEmergence.csv`: Emergence data from each replicate. Columns are:

**Site**: Numeric site distinguisher. Deprecated
**Site_Code**: Becomes Site_ID in code.
**Tray**: Numeric Tray (1-4)
**Tray_Code**: Becomes Tray_ID in code
**Class**: Land use class
**Month**: Numeric month (e.g. 8 = August)
**Day**: Day of month
*Exp_Day**: Day of experiment. Day 0 = day of hatching.
**Sex**: Sex of mosquito. M = male, F = female
**Num_Emerge**: Number of mosquitoes emerging from that tray on that day.

`AugustWingLength.csv`|`OctoberWingLength.csv` : wing lengths of female uninfected mosquitoes from each replicate.

**TrayCode**: Same as Tray_ID
**Date**: Date of emergence, written as MM/DD/YY
**Individual**: Number that corresponds to individual per tray and day
**Bars**: measurement in bars
**Conversion(bars/mm)**: To convert bars to mm

### clean

`emergenceTray.RData`: Mean development rate (1/Exp_Day) per tray.

`growthRates.RData`: Calculated growth rate by tray, via `03_model-calculations.Rmd`. Includes microclimate data.

`individuals.RData`: Individual level emergence for every potential female mosquito (50/tray). Event = 1 is emerged, Event = 0 is not emerged. Exp_Day is day of emergence.

`survivalTray.RData`: Number of female mosquitoes surviving until emergence per tray (numFSurv), out of 50. Not currently used directly in analysis.

`wingLength.RData`: Wing length of emerged, uninfected female mosquitoes. Includes day of emergence and associated temperature data (temperature from Exp_day=1 until day of emergence).

## infections

### raw

`augustDengue.csv`: Individual level infection data from summer replicate, 2016.
`octoberDengue.csv`: Individual level infection data from fall replicate, 2016.

Columns are:

**Individual**: Code for mosquito
**DPI**: Days post infection that infection was assessed
**Body, Head, Saliva**: Infection status at each level (Y=positive, N=negative, C = contaminated)
**WingLength**: wing length of mosquitoes in 'bars'
**conversion (mm/bar)**: conversion from bars to mm

### clean

`efficiencyPlot.RData`: Infection data summarized as effiency data. Used to plot.

`infectionSummary.RData`: Infection data summarized by land class and block, including mean and SE. Used primarily to plot.

`infectionSummaryLong.RData`: Infection data summarized by land class and block, including mean and SE. Used primarily to plot. In long format.

`seasonInfection.RData`: Cleaned individual level infection data for 21 dpi only. Covers both replicates and includes wing length of infected mosquitoes.

`VecCapacity.RData`: Calculated vectorial capacity. Created via `03_model-calculations.Rd`.

## microclimate

### raw

This contains all of the microclimate data as downloaded from data loggers.  Columns are as follows:

**Date**: Date and Time of recording
**Temp**: Temperature (C)
**RH**: Relative Humidity (%)
**Pot_ID**: Data Logger ID
**Site_ID**: ID of Site, R = rural, S = suburban, U = urban
**Class**: Class

### clean

`2016TrailsAdult.csv`: Cleaned climate data, matched to the data logger. Malfunctioned data loggers have been dropped and data is still in 10 min frequency. Also as `RData` file.

`climateTraySummary.RData`: Microclimate data by tray, averaged over the whole study period by season.

`fallInfectionClimate.RData`: Microclimate data by site, weighted and averaged to match the infection data. Includes mean and standard error.

`summerInfectionClimate.RData`: Same as fall, but for summer replicate.

### Other

`2015Climate.csv`: Climate data from 2015 experiment for comparison
`trayLoggerID.csv`: Matches the data logger ID to the tray ID in experiment

## spatial

These files are used primarily to create the map of sites.

**sites.csv**: lat-long coordinates of sites

**impSurfFocal.tif**: map of impervious surface used to categorize sites

**accMap.png**: outline of Athens Clark county

**cb_2016_us_county_500k.shp**: counties of the US