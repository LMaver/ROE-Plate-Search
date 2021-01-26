import numpy as np
from numpy import sin, cos
import sys
import math
from astropy.time import Time
from astropy.table import QTable
from astroquery.jplhorizons import Horizons
import time

class plate(object):
    def __init__(self, prefix, number, suffix,
                 survey_code,non_survey_code, ID,
                 RA, DEC, date, exposure_start, emulsion,
                 filter, exposure_duration, grade, prism_code):
        self.prefix = str(prefix)
        self.number = int(number)
        self.suffix = str(suffix)
        self.s_code = str(survey_code)
        self.non_s_code = str(non_survey_code)
        self.ID = str(ID)
        
        self.RA = float(RA)
        self.DEC = float(DEC)
        self.date = np.array([int(date[:2]), int(date[2:4]), int(date[4:])])
        
        if exposure_start[:-2] == '  ':
            hrs_int = 0
        else:
            hrs_int = int(exposure_start[:-2])
        
        self.start_time = np.array([hrs_int, int (exposure_start[-2:])]) #[hh, mm]
        self.emulsion = str(emulsion).replace(" ", "")  
        self.filter = str(filter).replace(" ", "")
        self.duration = (int(exposure_duration)/10)/60 #Plate cat in mmmt (t:tenth of a minute) convert to hrs
        self.grade = str(grade)
        self.p_code = str(prism_code)
        
        
        if self.date[1] < 10:
            month_str = '0' + str(self.date[1])
        else:
            month_str = str(self.date[1])
        
        if self.date[2] < 10:
            day_str = '0' + str(self.date[2])
        else:
            day_str = str(self.date[2])
    
        
        self.date_str = str(yr_2_to_4(self.date[0])) + '-' + month_str + '-' + day_str
        #print(self.date_str)    
        mjd_0 = Time(self.date_str, format = 'iso').mjd
        
        
        
        mid_point_LST = self.start_time[0] + float(self.start_time[1])/60. + float(self.duration)/2.
        long = 149.07
        GST = mid_point_LST - (long/15.)
        UT = GST2UT(mjd_0, GST)
               
        self.MJD = mjd_0 + UT/24.
        

        #Estimates for the Limiting Magnitude of the Plate taken from the UKST Handbook
        """
        ---------------------!!!!!!!!!!!!!!!---------------------
        THIS IS NOT EXTENSIVE SOME PLATES WITH EMULSION 4415 
        WILL RETURN A LIMITING MAGNITUDE OF ZERO(OR LESS)
        AS THEIR FILTERS(OR LACK OF) HAVE NOT BEEN ACCOUNTED FOR
        ---------------------!!!!!!!!!!!!!!!---------------------
        """
        #
        self.limMag = 0.
        if (((self.emulsion == 'IIaO') or(self.emulsion == 'IIaO') or 
             (self.emulsion == 'IIIaJ') or (self.emulsion == '4415')) 
            and
            (self.filter =='NONE' or  self.filter == 'UG1')):
            
            Nom_ext_t = 180./60. #Nominal Exposure time in mins
            Nom_depth = 21.0    #Nominal Plate Depth
            
        elif ((self.emulsion == 'IIaO') or self.emulsion == 'IIaO' or 
              self.emulsion == 'IIIaJ' or  self.emulsion == 'IIaD'):
            Nom_ext_t = 1.
            if self.emulsion == 'IIIaJ':
                Nom_depth = 22.5
            else:
                Nom_depth = 21.0
        
        elif self.emulsion == 'IIIaF':
            if self.filter == 'OG590':
                Nom_ext_t = 1.
            else:
                Nom_ext_t = 90./60.
                
            Nom_depth = 21.5
        
        elif self.emulsion == '4415':
            if self.filter == 'OG590':
                Nom_ext_t = 1.
                Nom_depth = 22.5
            else: 
                Nom_ext_t = 180./60.
                Nom_depth = 21.5
        
        elif (self.emulsion == 'IVN'):
            Nom_ext_t = 90./60.
            Nom_depth = 19.5
        else:
            Nom_ext_t = self.duration
            Nom_depth = 0.
        
        
        if self.duration != 0.:
            self.limMag = Nom_depth + 1.25*math.log(self.duration/Nom_ext_t)/math.log(10)
             
        else:
            self.limMag = 0.
        
        

    def __str__(self):
        return(str(self.number) + ' ' + str(self.RA) + ' ' 
               + str(self.DEC)  + ' ' + str(self.date) + ' ' 
               + str(self.limMag))


def plate_from_line(line):
    """
    Generates a plate object from a line of the machine readable catalogue
    :param file_handle: the handle of the plate catalogue
    """
    Data = str(line)
    Data = Data[:-1]
    #print(Data)
    prefix = Data[0:2]
    number = Data[2:7]
    suffix = Data[7]
    survey_code = Data[8:11]
    non_survey_code = Data[11:15]
    ID = Data[15:20]
    RA_hrs = float(Data[20:22])
    RA_mins = float(Data[22:24]) +(float(Data[24])/10.)
    RA = HMS_to_deg(RA_hrs, RA_mins)                    #Astropy methods to do this???
    DEC_deg = float(Data[25:28])
    #print(DEC_deg)
    DEC_min = float(Data[28:30])
    #print(DEC_min)
    DEC =  deg_min_sec_to_deg(DEC_deg, DEC_min)         #Probs this too, guess it doesn't matter
    #print(DEC)
    exposure_date = Data[30:36]
    exposure_start = Data[36:40]
    emulsion = Data[40:46]
    filter = Data[46:52]

    exposure_duration = Data[52:56]
    if len(Data) > 55:
        if suffix != 'P' :
            grade = Data[56:]
            prism_code = 'N/A'
        else:
            grade = Data[56:61]
            prism_code = Data[61:]
    else:
        grade = 'N/A'
        prism_code = 'N/A'

    return plate(prefix, number, suffix, survey_code, non_survey_code,
                 ID, RA, DEC, exposure_date, exposure_start, emulsion,
                 filter, exposure_duration, grade, prism_code)
          
def plate_list(filename):
    """
    Create a list of plates given the filename of the plate catalogue
    :param filename: String name of plate catalogue file
    :return: a list containing all the plates as plate objects
    """
    plate_list = []
    plate_catalogue = open(filename, 'r')

    for line in plate_catalogue:
        
        plate_list.append(plate_from_line(line))

    plate_catalogue.close()

    return plate_list

class ephemeris(object):
    def __init__(self, date, JD,  RA, DEC, mag, RA_unc, DEC_unc):
        '''
        :param date: in the form YYYY-MON-DD is converted to a [YY,MM,DD] array
        :param RA: in degrees
        :param DEC: in degrees
        :param mag: the apparent visual magnitude
        '''
        self.RA = float(RA)
        self.DEC = float(DEC)
        self.date = np.array([int(date[2:4]), int(month_to_num(date[5:8])), int(date[9:])])
        if self.date[1] < 10:
            month_str = '0' + str(self.date[1])
        else:
            month_str = str(self.date[1])
        
        if self.date[2] < 10:
            day_str = '0' + str(self.date[2])
        else:
            day_str = str(self.date[2])
    
        
        self.date_str = str(yr_2_to_4(self.date[0])) + '-' + month_str + '-' + day_str
        
        
        self.MJD = float(JD) - 2400000.5
        self.Vmag = float(mag)
        self.RA_unc = float(RA_unc)
        self.DEC_unc = float(DEC_unc)
        

    def __str__(self):
        return str(self.date) + ' ' + str(self.RA) + ' ' + str(self.DEC) + ' ' + str(self.Vmag)
    

def month_to_num(month):
    """
    Converts the three letter name of a month to its numerical value
    jan = 1
    dec = 12
    :param month:
    :return:
    """
    month_num = 0
    if month == 'Jan': month_num = 1
    elif month == 'Feb': month_num = 2
    elif month == 'Mar': month_num = 3
    elif month == 'Apr': month_num = 4
    elif month == 'May': month_num = 5
    elif month == 'Jun': month_num = 6
    elif month == 'Jul': month_num = 7
    elif month == 'Aug': month_num = 8
    elif month == 'Sep': month_num = 9
    elif month == 'Oct': month_num = 10
    elif month == 'Nov': month_num = 11
    elif month == 'Dec': month_num = 12
    else: print('error in date conversion')
    return month_num
  

def GST2UT(MJD, GST):
    """
    

    Parameters
    ----------
    MJD : the MJD at 0h on the greenwich calendar date
    GST : the GST in decimal hours

    Returns
    -------
    UT in hours

    """
    JD =  MJD + 2400000.5
    S = JD - 2451545.0
    T = S/36525.0
    T0 = 6.697374558 + (2400.051336 * T) + (0.000025862*(T**2))
    while T0 > 24 or T0 <0:
        if T0 > 24:
            T0 -= 24
        elif T0 < 0:
            T0 += 24
    
    A = GST - T0
    
    while A > 24 or A <0:
        if A > 24:
            A -= 24
        elif A < 0:
            A += 24
    
    UT = 0.9972695663*A
    
    return UT


def yr_2_to_4(year):
    """
    

    Parameters
    ----------
    year : int
        the year as 2 digits
        eg. 72 = 1972
            2 = 2002
    Returns
    -------
    the year as 4 digits

    """
    #Adjust year to include first 2 digits
    if year <= 2:
        year += 2000
    else:
        year += 1900
    
    return year
    
def HMS_to_deg(hours, minutes, seconds = 0.):
    """
    Converts hours mintutes seconds and converts them into degrees Target
    :param hours: 15 degrees
    :param minutes: 1/60 of an hour
    :param seconds: 1/60 of a second
    :return:
    """
    #print('hours %f, minutes %f, seconds %f'%(hours, minutes, seconds))
    return(15*hours + 15*minutes/60 + 15*seconds/(60**2))


def deg_min_sec_to_deg(degrees, minutes, seconds = 0.):
    """
    takes degrees minutes seconds and converts them into degrees
    :param degrees: int
    :param minutes: arcminutes 1/60 of a degree
    :param seconds: arcseconds 1/60 of an arcminute
    :return: the total of the three values in degrees
    """

    return(degrees +(minutes/60) + (seconds/(60**2)))


def s_to_tp(RA_deg, DEC_deg, RA_0_deg, DEC_0_deg):
    """
    transforms spherical cordinates to rectangular coordinates on the tangent 
    plane projection centred at tangent point RA_0 DEC_0
    
    Parameters
    ----------
    RA: float
        the Right Ascension of the target point (the ephemieris of the body)
        in degrees
    DEC: float
        the declination of the target point (the ephemeris) in degrees
    RA_0:float
        The Spherical RA coord of the tangent point(the plate centre)
        in degrees
    DEC_0:float
        the Spherical DEC coord of the tangent point(the plate centre)
        in degrees

    Returns 
    -------
    xi: float
        Horizontal Rectangular coordinate on the tangent plane
    eta: float
        Vertical Rectangular coordinate on the tangent plane
    j: int
        status: 0 = OK, star on tangent plane
                1 = error, star too far from axis
                2 = error, antistar on tangent plane
                3 = error, antistar too far from axis
    """
    #convert all values to radians
    RA = np.deg2rad(RA_deg)
    DEC = np.deg2rad(DEC_deg)
    RA_0 = np.deg2rad(RA_0_deg)
    DEC_0 = np.deg2rad(DEC_0_deg)

    RA_diff = RA - RA_0
    
   
    denom = sin(DEC)*sin(DEC_0) + cos(DEC)*cos(DEC_0)*cos(RA_diff)
    
    
    tiny = 1e-6
    
    if denom > tiny:
        j = 0
    elif denom >= 0.:
        j = 1
        denom = tiny
    elif denom > -1*tiny:
        j = 2
        denom = -1*tiny
    else:
        j = 3
        
         
    xi = (cos(DEC)*sin(RA_diff))/denom
    eta = (sin(DEC)*cos(DEC_0) - cos(DEC)*sin(DEC_0)*cos(RA_diff))/denom

    
    #Convert to degrees
    xi_deg = np.rad2deg(xi)
    eta_deg = np.rad2deg(eta)
    
    return xi_deg, eta_deg, j
   

def platesearch(body_name, plates, date_range, id_type = 'smallbody'):
    
    UKST = {'lon': 149.07,
            'lat': -31.27,
            'elevation': 1.13}
    
    
    obj = Horizons(id=body_name, id_type = id_type, location = UKST, epochs = date_range)

    ephem = obj.ephemerides(quantities = '1,9,36', refsystem = 'B1950') #Equinox B1950 same as plate Catalogue
    
    eph_list = []
    
    for row in ephem['datetime_str','datetime_jd', 'RA', 'DEC', 'V', 'RA_3sigma', 'DEC_3sigma']:
        eph_list.append(ephemeris(row[0][0:11], row[1], row[2], row[3], row[4], row[5], row[6]))
        

    found_plates = []
    
    for plt in plates:
        i = 0
        for eph in eph_list:
            if( ((eph.date[0]) == int(plt.date[0])) and
                ((eph.date[1]) == int(plt.date[1])) and
                ((eph.date[2]) == int(plt.date[2])) ):
                    
                    
                    eph_1 = eph_list[i+1]
                    scale = (plt.MJD - eph.MJD)
                 
                    RA = eph.RA + scale*(eph_1.RA - eph.RA)
                    DEC = eph.DEC + scale*(eph_1.DEC - eph.DEC)
                    
                   
                    
                    xi, eta, j = s_to_tp(RA, DEC, plt.RA, plt.DEC)
                    
                    if j == 0 and abs(eta) < 3.2 and abs(xi) < 3.2:
                        
                        x = ((xi*60.*60.) / 67.12) + 178  # Horizontal distance/coord from SW(SE?) corner in mm (67.14 arsec per mm)
                        y = ((eta*60.*60.) / 67.12) + 178 # Vertical distance/coord from SW corner in mm (67.14 arcsec per mm)
                        
                        max_RA_unc = max([eph.RA_unc, eph_1.RA_unc])
                        max_DEC_unc = max([eph.DEC_unc, eph_1.DEC_unc])
                        
                        
                        found_plates.append((plt.number, plt.grade, 
                                             plt.date, plt.MJD, 
                                             plt.RA, RA, 
                                             plt.DEC, DEC,
                                             max_RA_unc, max_DEC_unc, 
                                             plt.start_time, plt.duration, 
                                             plt.limMag, eph.Vmag, 
                                             x , y))
                        
            i += 1
    
    if found_plates == []:
        found_plates.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    
    return QTable(rows = found_plates, names = ('Plate #','date', 'MJD', 'Grade', 
                                                'Plate RA', 'Horizons RA', 'Plate DEC', 
                                                'Horizons DEC','RA Uncertainty', 
                                                'DEC Uncertainty', 'LST at Start', 
                                                'Exposure Duration', 'Limiting Mag', 
                                                'Target VMag',  'x', 'y' ))  
    #return found_plates
    
def search_plates(names, output_file, start_date = '1973-07-10', end_date = '2002-11-12'):
    """
    Given a list of names uses Plate reading to write a file containing plates
    found that may contain that body
    
    Parameters
    ----------    
    names : list
        a list of body names to be searched in the plate catalogue for

    survey : string
        The name of the survey to search through
        -UKST
        -ESO
        -POSS1
        -POSS2
    
    output_file: string
        The name the file to write search results to
    
    """
    
    output = open(output_file, 'w')
    output.write("Object Name, Plate Number, Plate Grade, Date, MJD, Plate RA, Horizons RA, Plate DEC, Horizons DEC, RA uncertainty, DEC uncertainty, Exposure Start Time, Exposure Duration, Limiting Magnitude, Object Magnitude, X-Coord, Y-Coord \n")
    filename =  'UKST_PlateCatalogue.txt'
    plates = plate_list(filename)
    date_range = {'start':start_date, 'stop': end_date, 'step': '1d'}

    
    n = 0
    interested_plates = 0
    maybes = 0
    un_calced = 0
    
    vector_platesearch = np.vectorize(platesearch, excluded=['plates'], otypes=[QTable])
    
    finds_tables = vector_platesearch(names, plates = plates, date_range =date_range)
    i = 0
    for table in finds_tables:
        for row in table:
            output.write(names[i])
            
            n += 1
            for item in row:
                output.write(', ' + str(item))
            
            if row[12] > row[13]:           #If LimMag > Object VMag
                interested_plates += 1 
            elif row[13] - row[12] <= 0.5:  #If VMag within 0.5 of LimMag
                maybes += 1
            elif row[12] < 10.:
                un_calced += 1              #Number of plates hit w/ undetermined LimMag
            
            output.write(' \n')
        
        print('Object search: ' + names[i])
        print('Plate Hits: ' + str(n))
        print('Plates with LimMag > Object VMag: ' + str(interested_plates))
        print('Plates with Object VMag within 0.5 of LimMag: ' + str(maybes))
        print('Plates with undetermined Limiting Magnitude: ' + str(un_calced))
        i+=1
    
    output.close()
    

        



def main():
    args = sys.argv
    
    #Check there's enough cmd line args
    if len(args) < 3:
        print("Usage 1: PlateSearch.py fromfile <input file> <output file>")
        print("Usage 2: PlateSearch.py <object names> <output file>")
        sys.exit(1)
        


  
    if args[1] == 'fromfile': #If no survey is specified and we're reading from file
        input_file = args[2]
        name_list =[]
    
        targets_file = open(input_file, 'r')
        for line in targets_file:
            name_list.append(str(line[:-1]))
        targets_file.close()
    
    else:
        name_list = []
        
        for i in range(1,len(args)-1):
            name_list.append(str(args[i]))
    output = args[len(args)-1]
    
    names = np.array(name_list)
    
    search_plates(names, output)
 
t_0 = time.time()
main()
print("Duration: %f"%(time.time()-t_0))