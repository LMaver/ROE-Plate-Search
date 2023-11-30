import numpy as np
from numpy import sin, cos
#import sys
import math
from astropy.time import Time
#from astropy.table import QTable
from astroquery.jplhorizons import Horizons
#import time
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
from pandastable import Table
#Load UKST to Dataframe.


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

def date2JD(date_string):
    #Method for converting YYYY-MM-DDTmm:MM:SS to Julian Dates
    #print(date_string)
    date, time = date_string.split(" ")
    D, M, Y  = date.split("/")
    hms =time.split(":")
    if len(hms) == 3:
        hour, minute, sec = hms
    elif len(hms) == 2:
        hour, minute, sec = hms[0], hms[1], 0.0 
    Y,M,D,hour, minute,sec  = int(Y),int(M),int(D),int(hour),int(minute),int(sec)
    #print(Y,M,D,hour, minute,sec)
    JDN = (1461 * (Y + 4800 + (M -14)/12))/4 + (367 * (M - 2 - 12 * ((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075
    
    JD = JDN + (hour-12)/24 + (minute/1440) + (sec/86400)

    MJD = JD - 2400000.5 #to convert to MJD
    #print(date_string, JD)
    return JD

datearray2JD = np.vectorize(date2JD)

def calc_lim_mag(emulsion, filter, duration):
    limMag = 0.
    
    if (((emulsion == 'IIaO') or(emulsion == 'IIaO') or 
        (emulsion == 'IIIaJ') or (emulsion == '4415')) 
        and (filter =='NONE' or  filter == 'UG1')):
        
        Nom_ext_t = 180./60. #Nominal Exposure time in mins
        Nom_depth = 21.0    #Nominal Plate Depth
        
    elif ((emulsion == 'IIaO') or emulsion == 'IIaO' or 
        emulsion == 'IIIaJ' or  emulsion == 'IIaD'):
        Nom_ext_t = 1.
        if emulsion == 'IIIaJ':
            Nom_depth = 22.5
        else:
            Nom_depth = 21.0
    
    elif emulsion == 'IIIaF':
        if filter == 'OG590':
            Nom_ext_t = 1.
        else:
            Nom_ext_t = 90./60.
            
        Nom_depth = 21.5
    
    elif emulsion == '4415':
        if filter == 'OG590':
            Nom_ext_t = 1.
            Nom_depth = 22.5
        else: 
            Nom_ext_t = 180./60.
            Nom_depth = 21.5
    
    elif (emulsion == 'IVN'):
        Nom_ext_t = 90./60.
        Nom_depth = 19.5
    else:
        Nom_ext_t = duration
        Nom_depth = 0.
    
    if duration != 0.:
        limMag = Nom_depth + 1.25*math.log(duration/Nom_ext_t)/math.log(10)
    else:
        limMag = 0.
    
    return limMag

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
    exposure_date_string = Data[30:36]
    date = [int(exposure_date_string[:2]), int(exposure_date_string[2:4]), int(exposure_date_string[4:])]
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

    if exposure_start[:-2] == '  ':
            hrs_int = 0
    else:
        hrs_int = int(exposure_start[:-2])
    
    start_time = np.array([hrs_int, int (exposure_start[-2:])]) #[hh, mm]
    duration = (int(exposure_duration)/10)/60 #Plate cat in mmmt (t:tenth of a minute) convert to hrs
    
    
    if date[1] < 10:
        month_str = '0' + str(date[1])
    else:
        month_str = str(date[1])
    
    if date[2] < 10:
        day_str = '0' + str(date[2])
    else:
        day_str = str(date[2])

    
    date_str = str(yr_2_to_4(date[0])) + '-' + month_str + '-' + day_str
    #print(self.date_str)    
    mjd_0 = Time(date_str, format = 'iso').mjd
    
    
    
    mid_point_LST = start_time[0] + float(start_time[1])/60. + float(duration)/2.
    long = 149.07
    GST = mid_point_LST - (long/15.)
    UT = GST2UT(mjd_0, GST)
            
    JD = 2400000.5 + mjd_0 + UT/24.

    limMag = calc_lim_mag(str(emulsion).replace(" ", ""), str(filter).replace(" ", ""), duration)


    plate_dict = {"Plate Number": number, 
                  "RA": RA, "DEC": DEC,
                  "Date": date_str, "JD": JD,
                  "Plate Grade": grade, "Plate Lim": limMag
                  }
    return plate_dict

def s2rec(RA_ephem, DEC_ephem, RA_obs, DEC_obs):
    """
    transforms spherical cordinates to rectangular coordinates on the tangent 
    plane projection centred at tangent point RA_0 DEC_0
    
    Parameters
    ----------
    RA_ephem: float
        the Right Ascension of the target point (the ephemieris of the body)
        in degrees
    DEC_ephem: float
        the declination of the target point (the ephemeris) in degrees
    RA_obs:float
        The Spherical RA coord of the tangent point(the plate centre)
        in degrees
    DEC_obs:float
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
    RA = np.deg2rad(RA_ephem)
    DEC = np.deg2rad(DEC_ephem)
    RA_0 = np.deg2rad(RA_obs)
    DEC_0 = np.deg2rad(DEC_obs)

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

def load_UKST():
    UKST_filename = "UKST_PlateCatalogue.txt"
    dict_list = []
    with open(UKST_filename, 'r') as cat:
        for line in cat:
             dict_list.append(plate_from_line(line))

    UKST_cat = pd.DataFrame(dict_list)
    
    return UKST_cat

def get_ephems(body, epochs):
        #method to retreive the ephemera of a given body in a given date range [epochs]
        #Coords of the Javalambre Observatory
        UKST = {'lon': 149.07,
            'lat': -31.27,
            'elevation': 1.13}
        
        obj = Horizons(id=body, location = UKST, epochs = epochs)
        ephem = obj.ephemerides(quantities = '1,9,36', refsystem = 'B1950') #Equinox B1950 same as plate Catalogue

        return ephem #ephem is a Astroquery pd.Datframe-esque object

def check_obs(JD, RA, DEC, ephems):
    #Check if RA & DEC are close to ephems on given date
    #print(JD, RA, DEC)
    #print(ephems)
    #Get the .5 JD prior to & post the observation
    JD_before = np.round(JD) - 0.5
    JD_after = JD_before + 1
    #print(ephems['RA'])#[ephems['datetime_jd']==float(JD_before)])

    #Pull out the RA,DEC from the closest ephemera prior to obs
    RA_before = ephems["RA"][ephems['datetime_jd']==JD_before].iloc[0]
    DEC_before = ephems["DEC"][ephems['datetime_jd']==JD_before].iloc[0]

    #Pull out the RA,DEC from the closest ephemeris after the obs
    RA_after = ephems['RA'][ephems['datetime_jd']==JD_after].iloc[0]
    DEC_after = ephems["DEC"][ephems['datetime_jd']==JD_after].iloc[0]

    #interpolate between the before & after to get approx RA,DEC at time of obs 
    scale = (JD - JD_before)
    RA_scaled = RA_before + scale*(RA_after - RA_before)
    DEC_scaled = DEC_before + scale*(DEC_after - DEC_before)
    
    #Calculate obj distance from centre
    #Trasform Coords, unsure if this is needed i can't justify why not, a plate and a digital image are identical for all intents
    H_dist, V_dist, j = s2rec(RA_scaled, DEC_scaled, RA, DEC)

    #H_dist = abs(RA - RA_scaled)
    #V_dist = abs(DEC - DEC_scaled)
    dist  = np.sqrt(H_dist**2 +  V_dist**2)
    #If the Ephemeris PosN is within 3 degrees of the Obs PosN return it as a match --!!-- Double check the field of view 
    if j==0 and (abs(H_dist)<3.2) and (abs(V_dist) < 3.2):
        x = ((H_dist*60.*60.) / 67.12) + 178  # Horizontal distance/coord from SW(SE?) corner in mm (67.14 arsec per mm)
        y = ((V_dist*60.*60.) / 67.12) + 178
        dist = math.sqrt(x**2 + y**2) 
        return True, JD, dist, x, y
    else: return False, JD, 0.0, 0.0, 0.0

#Vectorize the method
check_array = np.vectorize(check_obs, excluded=["ephems"])

def body_search(body, start_date = '1973-07-10', end_date = '2002-11-12'):
    date_range = {'start':start_date, 'stop': end_date, 'step': '1d'}
    
    ephems = get_ephems(body, date_range)['datetime_jd', 'RA','DEC', 'V' ].to_pandas()
    results = check_array(UKST_df['JD'], UKST_df['RA'], UKST_df['DEC'],  ephems = ephems)
    
    results_df = pd.DataFrame({'Match': results[0],"JD": results[1], 
                            "Dist From Centre": results[2],"x": results[3], "y": results[4]})
    matches_df = results_df[['JD', 'Dist From Centre', 'x', 'y']][results_df['Match']]
    avg_Vmag = ephems['V'].mean()
    final = UKST_df.merge(matches_df, how = 'right', on = 'JD')
    final["Good Match"] = final["Plate Lim"] > avg_Vmag + 1.5
    N_goodMatches = len(final[final["Good Match"]].index)
    N_unsureMatches = len(final[final["Plate Lim"] == 0].index)
    return final, N_goodMatches, N_unsureMatches, avg_Vmag

def run_window():
    def multi_search(body_list):
        N_of_matches = []
        N_goodMatches = []
        N_unsureMatches = []
        Vmags  =[]
        i = 1
        for body in body_list:
            lbl_progress.config(text= "Searching for: %s (%s/%s)"%(body, i, len(body_list)))
            summary_df = pd.DataFrame({"Body Name": body_list[:i-1], "# of Plate Matches": N_of_matches,"# of Good Matches": N_goodMatches,"# w/o Plate Lim": N_unsureMatches, "Avg VMag of Body": Vmags})
            summary_tb.model.df = summary_df
            stb.redraw()
            window.update()
            
            body_matches, N_good, N_unsure, avg_Vmag = body_search(body)
            N_of_matches.append(len(body_matches.index))
            N_goodMatches.append(N_good)
            N_unsureMatches.append(N_unsure)
            Vmags.append(avg_Vmag)
            body_matches.to_csv(body+"_platematches.csv")
            i+=1

        summary_df = pd.DataFrame({"Body Name": body_list[:i-1], "# of Plate Matches": N_of_matches,"# of Good Matches": N_goodMatches,"# w/o Plate Lim": N_unsureMatches, "Avg VMag of Body": Vmags})
        summary_tb.model.df = summary_df
        stb.redraw()

        df = pd.read_csv(body_list[0] + "_platematches.csv", index_col=0)            
        df_tb.model.df = df
        pt.redraw()

        lbl_progress.config(text= "Completed Search for %s bodies"%(len(body_list)))

    def press_multi_search(var = None):
        #Run multiple body search
        #Display Summary Results
        body_list = multi_entry.get().split(",")
        if len(body_list) != 0:
            multi_search(body_list)
            
    def select_file(var = None):
        # On browse btn press
        # Select file from folder
        # Run multiple bodies search
        # Display Summary Search
        # Create Dropdown Selection?
        filetypes = (
            ('txt files', '*.txt'),
            ('All files', '*.*')
            )

        filename = fd.askopenfilename(
            title='Open files',
            initialdir='./',
            filetypes=filetypes)
        
        body_list = []
        with open(filename, 'r') as f:
            for line in f: body_list.append(line.replace("\n", ""))
        
        if len(body_list) > 0: multi_search(body_list)         

    #Toplevel
    window = tk.Tk()
    window.title("Plate Search")
    window.resizable(width=True, height=True)
    
    frm_display = tk.Frame(master = window, height = 200,width = 1000)
    left_frm = tk.Frame(window)
    summary_frm = tk.Frame(window, height = 20)
    
    frm_display.grid(row=2, column = 0, padx=10)
    left_frm.grid(row = 0, column = 0,  padx=10)
    summary_frm.grid(row=1,column=0,  padx=10)
    
    #Display table
    df = pd.DataFrame()
    df_tb = pt = Table(frm_display, dataframe=df, showtoolbar=True, showstatusbar=True, width = 700)
    df_tb.grid(row = 0,column = 0)    
    pt.show()

    
    #Multi-Body Search
    multi_lbl = tk.Label(left_frm, text = "For Multiple bodies enter names comma separated or select a file")
    multi_entry = tk.Entry(left_frm, width=40)
    
    #brs_btn = tk.Button(left_frm, width = 10, text="Select File", command = select_file)
    multi_lbl.grid(row = 0, column  = 0, columnspan = 2)
    multi_entry.grid(row = 1, column = 0)

    lbl_progress = tk.Label(left_frm, text = "")
    lbl_progress.grid(row=2, column = 0, sticky = 'w', columnspan = 2)

    #Button
    btn_multi_go = tk.Button(
    master = left_frm,
    width = 10,
    text = "Select File",
    command = lambda: select_file(None))
    #brs_btn.grid(row = 1, column = 1)
    btn_multi_go.grid(row = 1, column = 1, sticky = "w")
    multi_entry.bind('<Return>', press_multi_search)
    
    sdf = pd.DataFrame()
    summary_tb = stb = Table(summary_frm, dataframe = sdf, 
                             showtoolbar=False, showstatusbar=False, 
                             height = 100, width = 720)
    summary_tb.grid(row=0, column = 0)
    stb.show()


    window.mainloop()

def main():
    #Load Catalougue
    global UKST_df
    UKST_df = load_UKST()
    run_window()

main()