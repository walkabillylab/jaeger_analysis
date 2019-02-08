## Authors
* [Javad Rahimipour Anaraki](https://github.com/jracp)
* [Daniel Fuller](https://github.com/walkabillylab)


## Introduction

Code to clean data from [Jaeger Oxycon Pro](https://www.ntnu.edu/documents/221360533/221362168/oxypro.pdf/91b256f4-97d7-4408-8eec-53869132a03d) metabolic cart based on the paper "[Methods Used to Process Data from Indirect Calorimetry and Their Application to VO2MAX](https://www.asep.org/asep/asep/Robergs3.pdf)" Robert A Robergs and Angus F Burnett. Journal of Exercise Physiology. 6(2), 2003. 


## Preliminaries

Place the code next to the participants folder which contain a file called **timing.csv** and a CSV data file generated by Jaeger. 

* data
    + Jaeger
        + 301
            + timing.csv
            + JaegerDataFile.csv
        + 302
            + timing.csv
            + JaegerDataFile.csv
        + ...
    + intervals.csv

The **timing.csv** should contain the following columns:
 
* `Task name` stores the task done by participant
* `Start time` is start date and time (in ####-##-## ##:##:## format)
* `End time` is end date and time (in ####-##-## ##:##:## format)

The **data** folder also contains a file called **intervals.csv** which stores information about each participants as follows:
 
* `herox` stores participants alias names
* `kit` is the package number
* `phone` is phone id
* `watch` is Apple watch id
* `fitbit` is Fitbit id
* `geneactiv` is GENEActiv id
* `snesedoc` is SenseDoc id
* `start` is start date and time (in ####-##-## ##:##:## format)
* `end` is end date and time (in ####-##-## ##:##:## format)
* `userid` is a unique id for each participant
* `wrist` is 1 if Apple Watch and GENEActiv are on the same wrist, otherwise 0
* `age`, `gender`, `weight` and `height` are demographic data
* `street`, `city`, `postal` are address information for each participant.

**Note**: This file should be kept updated throughout the experiment.


## Run
Steps to run the code as in the follows:

### 1. Modifying required variables

* Modify `path` and `intrPath`, the Jaeger participants folder and **intervals.csv**, respectively
* Change the `uid` (e.g. 301) for choose which participant data should be processed
* Modify the `timeZone` variable to your local timezone

**Info**: More information on setting the `timezone` can be found [here](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). 

### 2. Run the code
Select all the code and run it using Ctrl + Enter.

### 3. Results
Based on the interpolation and spline cure fitting six plots are generated. For VO2, VCO2 and RER the generated curve and Bland-Altman figures are plotted and saved in each participant folder. Also, a CSV file to check the validity of the collected Jaeger data is generated and placed next to the plots.

