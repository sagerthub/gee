/*
Last update: 2022-10-06
Author: MRMPO (gis@mrcog-nm.gov), adapted from:
	Sean McCartney (sean.mccartney@nasa.gov)
	ARSET Training: Satellite Remote Sensing for Measuring Urban Heat Islands and Constructing Heat Vulnerability Indices
	August 2, 2022 - August 11, 2022

This code is free and open. 
By using this code you agree to cite the following reference in any publications derived from them:
NASA Applied Remote Sensing Training (ARSET) program

This example shows how to analyze and visualize Landsat surface temperature (ST) time series 
from Landsat 8 central New Mexico (USA) from a defined area of interest (aoi).

Parameters:
In: DATE_RANGE
    YEAR_RANGE
    STUDYBOUNDS
    DISPLAY
    aoi: delineated rectangle for area of interest
    
Out: clipped image of mean daytime surface temperature (ST).
*/


//****************** USER INPUTS *********************//
/*
Note:
To modify export, see lines 188-196.
To modify map window image, see lines 190-198.
*/
/*
Enter date as the day of the year.
e.g., July 1 (non-leap-year) is day 182 and Aug 31 is day 243.
*/
var start_day = 182;
var end_day = 243;

// Enter start and end years of range (inclusive).
var start_year = 2013;
var end_year = 2022;

/* DEFINE THE STUDY AREA 

The following coordinates define a square around the MRCOG area,
from NW of Sandoval County to SE of Torrance County.
Alternative methods below. Note that in testing,
other shapes failed to export.
*/
var aoi = 
    ee.Geometry.Polygon(
        [[[-107.655533, 36.242138],
          [-105.265001, 36.242138],
          [-105.265001, 34.233749],
          [-107.655533, 34.233749]]], null, false);

// Alternative AOI Methods
// Use the "Draw a rectangle" tool in the map window to define your 
// area of interest (aoi) for filtering the image collection.
// Or use an Asset (e.g., imported shapefile).
// Under Imports, rename variable to aoi. Or
// change aoi above or STUDYBOUNDS below to match import name.

var STUDYBOUNDS = aoi;

// Set the basemap to display as satellite.
Map.setOptions('SATELLITE');
/* Other basemap options are below
Map.setOptions('HYBRID');
Map.setOptions('ROADMAP');
Map.setOptions('TERRAIN');
*/

// Center the map view at defined coordinates (longitude/latitude) with the given zoom level.
// Zoom level ranges from 0 to 24, with 0 being global and 24 being the smallest region possible.
// Latitude must be within [-85, 85].
// Adjust the long/lat for your own aoi.
/* Option 1: Center at Big I in Albquerque, zoom level 9.
	Comment (//) one option so only the other is used. */
// Map.setCenter(-106.629260, 35.105125, 9);
/* Option 2: Center at centroid of AoI.
	Comment (//) one option so only the other is used. */
var geom = STUDYBOUNDS.centroid(0.001);
var coords = geom.coordinates();
var firstC = coords.get(0);
var lastC = coords.get(-1);
var xCoord = firstC.getInfo();
var yCoord = lastC.getInfo();
Map.setCenter(xCoord, yCoord, 9);
/*
Reference for Option 2:
https://gis.stackexchange.com/questions/382343/extracting-geometry-coordinate-from-a-feature
https://stackoverflow.com/questions/53629403/google-earth-engine-ee-number-to-integer
*/


//****** DO NOT MODIFY ******//

var DATE_RANGE = ee.Filter.dayOfYear(start_day, end_day);
// Assign a variable to filter years from 2010 – 2022.
// Adjust the YEAR_RANGE for your own UHI study.
var YEAR_RANGE = ee.Filter.calendarRange(start_year, end_year,'year');

// Assign a variable to display images in the map window
var DISPLAY = true;

// Assign a variable to the sensor-specific bands unique to each Landsat mission.
var LC08_bands = ['ST_B10', 'QA_PIXEL']; // Landsat 8 surface temperature (ST) & QA_Pixel bands

//****************** CLOUD MASK FUNCTION *****************//
// Create a function to mask clouds and cloud shadows based on the QA_PIXEL band of Landsat 8 & 9
// For information on bit values for the QA_PIXEL band refer to: 
// https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2#bands
function cloudMask(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3)
    .or(qa.bitwiseAnd(1 << 4));
  return image.updateMask(mask.not());
}

/* Assign variables to import the Landsat Collection 2, Tier 1, Level 2 image collections, selecting 
the ST and QA_PIXEL bands, and spatially filtering the image collection by your aoi. */
var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .select('ST_B10', 'QA_PIXEL')
  .filterBounds(STUDYBOUNDS)
  .filter(DATE_RANGE)
  .filter(YEAR_RANGE)
  .map(cloudMask);

// Filter the collections by the CLOUD_COVER property so each image contains less than 20% cloud cover.	
var filtered_L8 = L8.filter(ee.Filter.lt('CLOUD_COVER', 20));

// Use print statements to print the argument to the console.
print(filtered_L8, 'Landsat 8 ST');

/* Create a funtion using Landsat scale factors for deriving ST in Kelvin and Fahrenheit.
For more information on ST scale factors, refer to:
https://www.usgs.gov/landsat-missions/landsat-collection-2-level-2-science-products */
function applyScaleFactors(image) {
  var thermalBands = image.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).multiply(1.8).add(32);
  return image.addBands(thermalBands, null, true);
}

// Use print statements to print the argument to the console.
print(filtered_L8, 'Landsat ST (Fahrenheit)');

// Define a variable to apply scale factors to the filtered image collection.
var landsatST = filtered_L8.map(applyScaleFactors);

// Use a print statement for tracking your progress in the console tab.
print("... Computing mean ST across image collection");

//****************** CALCULATE MEAN SURFACE TEMPERATURE *****************//
// Define a variable to calculate mean ST for each pixel geography 
// throughout the filtered image collection.
var mean_LandsatST = landsatST.mean();

// Extract just ST band for export (i.e., remove QA band).
var one_band = mean_LandsatST.select('ST_B10');

// Define a variable to use the clip funtion to subset your imagery to the aoi.
var clip_mean_ST = one_band.clip(STUDYBOUNDS);

// Use a print statement to print the argument to the console.
print(clip_mean_ST, 'Mean ST clipped to study area');

// Define a variable to select the temperature band.
var values_ST = clip_mean_ST.select("ST_B10"); 

// Define a variable to output a histogram of mean ST values within your aoi.
// Histogram is useful for seeing min/max values, which can be used for visualization (see below).
var histogram_ST_values = ui.Chart.image.histogram({
  image: values_ST, 
  region: STUDYBOUNDS, 
  scale: 30, 
  maxPixels: 1e9
});

// Use  print statement to output the histogram values of mean ST to the console tab.
// Histogram is useful for seeing min/max values, which can be used for visualization (see below).
 print(histogram_ST_values);

// Add the image to the map window, defining min/max values, a palette for 
// symbology, assign a name to the visualization, and display the result.
// Uses 50 deg F and 140 deg F as the limits for color classification.
// Those values may be changed if needed.
Map.addLayer(clip_mean_ST, {
  bands: "ST_B10", 
  min: 50, max: 140, 
  palette: ['blue','white','red']}, "ST", DISPLAY);

Export.image.toDrive({
  image: clip_mean_ST,
  description: 'MeanST_MRCOG_2013_2022_JulAug',
  folder: 'SUHI',
  fileNamePrefix: 'MeanST_MRCOG_2013_2022_JulAug',
  scale: 30,
  maxPixels: 1e12,
  fileFormat: 'GeoTIFF'
});


//****************** CALCULATE SUHI *****************//
/* To calculate the surface urban heat island (SUHI) intensity for your aoi...

Subtract the mean urban temperature from the mean rural temperature.
(ST_B10_mean [Urban]) - (ST_B10_mean [Rural]) = SUHI intensity

The intensity of the heat island is the simplest quantitative indicator of the 
thermal modification imposed by a city upon the geography in which it is situated, and 
of its relative warming in relation to the surrounding rural environment. */