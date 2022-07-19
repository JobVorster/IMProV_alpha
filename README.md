# IMProV_alpha
(In Development) Alpha version of Interactive Maser Proper motion for VLBI. This code aims to be a general toolkit for high resolution maser astrometry, after calibration and imaging. 

Note that this is a pre-release version of the code. It has been used to calculate proper motions from maser positions for multiple epochs, but does not yet have the full functionality planned for the project.

Functionality of this release:  <br /><br />
  script: "one_feature_identification.py":  <br />
    Interactive plots of maser positions, and spectra. <br />
    Single Gaussian fitting on the spectral profile of selected regions.<br />
    Saving of Gaussian fitting (selection region, and Gaussian fitting parameters, with uncertainties on the fits).<br /><br />
  script: "two_gaussian_to_feature_maps.py":<br />
    Calculates the flux weighted centroid of the maser features identified in the previous script.<br />
    Calculates the uncertainties with a monte-carlo simulation. <br /><br />
  script: "three_calc_pm.py":<br />
    Identify persistent maser features and calculate their relative proper motions with linear regression.<br />
    Save the proper motion calculations.<br />
    <br /><br />

How to use the code:<br />

script: "one_feature_identification.py":<br />
  Run the program on the terminal as the user runs python scripts on their machine.<br />
  Copy and paste the link into your browser.<br />
  
  The user can upload a .csv file with the correct column names (outlined below), and the maser spots in the csv files will be plotted. The user can then pan, zoom and zoom with the tools at the top left corner of the left plot. The user can also use "box select" shown with a rectangle in the top left to select regions. The spectra of the selected region will then be plotted on the central plot. The user can then click on the "Do Gaussian Fit" checkbox to fit one Gaussian to the spectral profile. The user can then press the "SAVE CURRENT GAUSS FIT" button to save the region and Gaussian fit parameters. It will be shown in the table on the right hand side. The user can also press the "Export" button to save the Gaussian fitting and regions for this epoch.

This script does not need any modification to work (other than the dependancies), except that the file you upload must:<br />
  Be ONE file.<br />
  Ending in ".csv" and be comma separated (no spaces).<br />
  With the following columns: "RA,RAERR,DEC,DECERR,VLSR,FLUX,DFLUX"<br />
    RA: Right ascension offset (in arcsec).<br />
    RAERR: Uncertainty on right ascension (in arcsec).<br />
    DEC: Declination offset (in arcsec).<br />
    DECERR: Uncertainty on declination (in arcsec).<br />
    FLUX: Peak intensity of maser spot (in Jy/beam).<br />
    DFLUX: Uncertainty in peak intensity of maser spot (in Jy/beam).<br />
 If the user has made a mistake or wants to process a new dataset, they have to restart the program (at this stage you can not refresh the saved Gaussian fits). <br />
 <br />
 script: "two_gaussian_to_feature_maps.py":<br />
  Before the user runs this script. They first need to add two things (it is shown clearly:<br />
    The path to the folder containing the original data tables (with the same columns as outlined above).<br />
    The path to the folder containing the Gaussian fitting parameters from the previous script.<br />
   The user can then run the script as they run python scripts on their machine.<br />
  The script will output csv files that have the appendix "feature" in the folder containing the maser spot tables.<br />
  
  script: "three_calc_pm.py":<br />
  Before the user runs this script, they need to add two quantities in the script (it is shown clearly at the top of the script).<br />
    The source declination (in decimal degrees) e.g. -35.54461<br />
    A list of the observation dates (in decimal years) e.g. [2014.551, 2014.669, 2015.336]. The dates should correspond to the dates of the maser spot tables in question.<br />
    After the user runs the script. They should copy and paste the link into their browser.<br />
    They can then upload MULTIPLE files. These files should be the output of the second script, outlined above.<br />
    On the left side, the plot will show the positions of the maser features' intensity weighted centroid. The color indicates the epoch of the observation.<br />
    The user can mouse over the individual data points.<br />
    The user can then select persistent maser features (ONLY ONE FROM EACH EPOCH) by clicking on the first data point, and then holding "SHIFT" to select more data points. <br />
    When the user has selected the data points of a single maser feature over multiple epochs the user can press the "SAVE FEATURES" button. It will calculate and show the proper motion for that feature. A red arrow will show the proper motion of that feature, and the maser feature positions will be hidden.<br />
    The "CALCULATE PROPER MOTIONS" button does nothing at this stage.<br />
    Below the buttons, on the right hand side. Is the value of the "Likeness Parameter". This is an experimental number that calculates how much alike the selected maser features are based on their spectral profiles. If you have any questions about this parameter, please contact me personally. You can ignore this parameter if you like, it does not affect the proper motions.<br />
    After you have selected the persistent features, and calculated the proper motions, you can "EXPORT" the table with the proper motions. It contains the proper motions in arcsec/year. <br />
    
    
At this stage, the program only calculates linear proper motions, and does not calculate parallax, peculiar motion and galactic rotation.
These are being worked on. Further, the program cannot fit more than one Gaussian to the spectral profile. <br />

If you use this program, if anything is unclear in the documentation, or if you find any bugs, please send me an email at:<br />
jobvorster8@gmail.com<br />
This is a "pre-release" so I do expect many corrections and improvements will be needed.<br />
