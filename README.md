# ROE Plate Catalogue Search

A tool for searching the catalogue of UKST plates located at
the Royal Observatory of Edinburgh for small bodies

# Update : Nov '23
'ImprovedUKSTsearch.py' runs faster and has a GUI for ease of use. Run Program, enter small body name(s) and table of plate match results are produced. 

# Included Files
- PlateSearch.py 
- UKST_Catalogue.txt

# Usage
Run from the command line

1. Searching for a single body
		python PlateSearch.py <Object Name> <Output File>
		
		eg. To search for Pluto and output the results to a file named 'pluto_results.txt', run:
		
		python PlateSearch.py Pluto pluto_results.txt
		
2. Searching for multiple bodies
	There are two ways to search for multiple bodies in a single search
	
	- Directly in the Command line
		python PlateSearch.py <Object 1> <Object 2> ... <Object N> <Output>
		
		eg. To Search for both Pluto and Eris and output to a file named 'dwarf_planets.txt'
		
		python PlateSearch.py Pluto Eris dwarf_planets.txt
		
	- From a text file
		python PlateSearch.py fromfile <Input File> <Output File>
		
		where < Input File > is the name of a file where each line is the name of a body to be searched for
		
		eg. to search for bodies named in a file called 'target_bodies.txt' and output to a file named 'plate_results.txt':
		
		python PlateSearch.py fromfile target_bodies.txt plate_results.txt

It typically takes around 4 minutes to search the entire catalouge from 1973 to 2002 for a single body
There's certainly room for optimisation in the code and the date range to be searched can be adjusted in 
the search_plates module of PlateSearch.py if not the entire date range needs to be searched
		
# Results
Firstly for each object searched for the command line displays
	- the number of plate hits, 
	- the number of plates found where the object's visual magnitude is less than the limiting magnitude of the plate
    - the number of plates found where the object is no more than 0.5 greater than the limMag of the plate
	- the number of plates found where the limMag of the plate hasn't been calculated 

The results of the search are also written to a text file
Each line is a plate on which the desired object might be found on
The line contains:
	- The Name of the Target Object
	- The Plate Number
	- Plate Quality Grade
	- Date of the Plate Exposure
	- MJD of Plate Exposure
	- RA of the Plate
	- RA of the Object from JPL Horizons
	- DEC of the Plate
	- DEC of the Object from JPL Horizons
	- Uncertainty in the Horizons RA
	- Uncertainty in the Horizons DEC
	- The LST start time of the exposure as [HH, MM]
	- The duration of the exposure in decimal hours
	- The Limiting Magnitude of the Plate
	- The Visual Magnitude of the Object from JPL Horizons
	- The X-Coordinate
	- The Y-Coordinate

Delimited by a comma
	
The (X,Y) coord is the position of the object on the plate in mm, measured from the bottom-left corner of the back of the plate

Further details of the plate catalogue can be found at: https://www.roe.ac.uk/ifa/wfau/ukstu/pltcat.html
as well as details of the UKST and the plates in general at: https://www.roe.ac.uk/ifa/wfau/ukstu/telescope.html

# Extras
The file containing the machine readable plate catalogue is included (UKST_Catalogue.txt)
but can also be found through the above roe pltcat link
and directly at: https://www.roe.ac.uk/ifa/wfau/ukstu/ukst_catalogue.lis

Additionally, the POSS and ESO catalogues use the same format and as such can be searched using PlateSearch.py
The catalogues for these surveys can be found here: https://www.roe.ac.uk/~nch/ap4.html
and simply change the catalogue file name in search_plates module of PlateSearch.py
