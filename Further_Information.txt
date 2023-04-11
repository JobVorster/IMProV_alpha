Interactive Proper Motions for VLBI (IMProV)

IMProV is a set of three Python scripts that allow for the interactive identification and calculation of proper motions for maser features on single-epoch VLBI spot maps.
IMPRoV 0.5.0 - One Feature Identification

This is the first script of the program. Its aim is to interactively identify maser features on single epoch spot maps. You can upload a CSV file with the required columns ("['RA','RAERR','DEC','DECERR','VLSR','FLUX','DFLUX']") and use the Zoom and Select tools to navigate the data. Use the rectangle tool to select maser spots, and if they show a Gaussian spectral distribution, you can press the "Save features" button to save a rectangle (x0,x1,y0,y1) containing the specific data points. Once you have selected all the maser features you want to select, you can press the "Export" button to save a CSV file containing the rectangles for each feature. Do this for each epoch which you have maser spot maps.

Author: Job Vorster
Date: April 11, 2023
Requirements: Python 3, dash, numpy, pandas, matplotlib, scipy, datetime
Usage: python one_feature_identification.py
Interaction with the program is via GUI as described above.
Contact: jobvorster8@gmail.com
IMPRoV 0.5.0 - Two Feature Statistics

This is the second script of the program. Its aim is to calculate feature statistics of detected maser spots on VLBI images. The program requires the user to input the foldernames containing the regions of interest (ROIs) as well as the foldername containing the original data. The statistics of the detected maser spots within the ROIs are then calculated, including their flux weighted centroid and a Gaussian fit is done to the data. If the Gaussian fit does not succeed, the parameters are replaced by 999.

Author: Job Vorster
Date: April 11, 2023
Requirements: Python 3, numpy, pandas, matplotlib, scipy, skspatial, glob
Usage: python two_feature_statistics.py SPOTMAPS_FOLDERNAME FEATURE_BOXES_FOLDERNAME
Contact: jobvorster8@gmail.com
IMPRoV 0.5.0 - Three Calculation of Proper Motions

This is the third and final script of the program. Its aim is to interactively calculate proper motions with the output of the second script in this program. It opens up an interactive plot that you can click on individual data points. Hold SHIFT to select more than one data point. Make sure to select only ONE DATA POINT per epoch, and obviously use your judgement on which is the correct features to select for proper motions. The VLSR, PEAK and FWHM are shown to help you judge. If you have selected the maser features and you want to calculate the proper motions, press the "SAVE FEATURE" button. It will be saved into a table, and an arrow indicating the proper motion will appear on the plot. Once you have selected all the features and calculated their proper motions, you can export the table to a CSV.

Author: Job Vorster
Date: April 11, 2023
Requirements: Python 3, numpy, pandas, dash, scipy, matplotlib, datetime
Usage: python three_calc_proper_motions.py SOURCE_PARS_FILENAME
Contact: jobvorster8@gmail.com
