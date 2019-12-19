import os
from datetime import datetime, timezone, timedelta
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pyproj as proj
import numpy as np

## Created by Semme J. Dijkstra 8/28/2019
    #  Added method: read_hypack_raw_file, Semme J. Dijkstra 10/4/2019

class Position:
    """A Class for handling Position Data"""

    def __init__(self):
       
        # The data attributes
        self.times = list()
        
        # The geodetic coordinates - these are curvilinear so do not put them
        # in vectors, as that is a linear concept
        self.latitudes = list()
        self.longitudes = list()
        self.heights = list()
        self.qualities = list()
        self.n_sats = list()
        self.hdops = list()
        self.separations = list()
        self.corr_ages = list()
        self.corr_stations = list()
        
        # The projected coordinates are (usually) in Cartesian coordinates (which are linear), store them in numpy arrays
        self.proj_pos=np.array([])
        self.data_path=str()
        self.metadata = dict()
        self.metadata["geodetic_units"] = "rad"
        self.metadata["height_units"] = "m"
        self.metadata["proj_units"] = "m"
        self.metadata["geoid_name"] = None
        self.metadata["ellipsoid_name"] = None
        self.metadata["height_relative_to"] = None
        self.metadata["time_basis"] = "UTC"
        self.metadata["proj_str"] = None


    def __str__(self): 
        txt  = "Geoid Used             : %s\n" % (self.metadata["location_name"])
        txt += "Ellipsoid Used         : %s\n" % (self.metadata["datum_type"])
        txt += "Reference Surface Name : %s\n" % (self.metadata["datum_name"])
        txt += "Observation Time Basis : %s\n" % (self.metadata["time_basis"])
        txt += "Observations Units     : %s\n" % (self.metadata["units"])

        if len(self.times):
            txt += "Start Time             : %s\n" % (min(self.times))
        else:
            txt += "No time data present\n"
        
        # If there are no latitudes, logically there are no longitudes either
        
        if len(self.latitudes):
            txt += "Minimum latitude       : %.2f\n" % (min(self.latitudes))           
            txt += "Maximum latitude       : %.2f\n" % (max(self.latitudes))       
            txt += "Minimum longitude      : %.2f\n" % (min(self.longitude))           
            txt += "Maximum longitude      : %.2f\n" % (max(self.longitude))       
        else:
            txt += "No geodetic coordinates are present\n"
            
        # Test heights separately, as we may just have horizontal positioning
            
#        if len(self.heights):
#             txt += "Minimum height       : %.2f%s\n" % (min(self.heights), self.metadata["height_units"]         
#             txt += "Minimum height       : %.2f%s\n" % (min(self.heights), self.metadata["height_units"]) 

        return txt


    def read_jhc_file(self, fullpath):

        # Set the reference ellipsoid to WGS84

        self.metadata["ellipsoid_name"] = "WGS84"
        self.metadata["geoid_name"] = "EGM08"
        self.metadata["height_relative_to"] = "geoid"
        
        # Check the File's existence
        if os.path.exists(fullpath):
            self.data_path = fullpath
            print('Opening GNSS data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)

        # Open, read and close the file
        gnss_file = open(fullpath)
        gnss_content = gnss_file.read()
        gnss_file.close
        
        times=list();

        # Tokenize the contents
        gnss_lines = gnss_content.splitlines()
        count = 0  # initialize the counter for the number of rows read
        for gnss_line in gnss_lines:
            observations = gnss_line.split()  # Tokenize the string
            time = datetime.fromtimestamp(
                float(observations[5]), timezone.utc)
            times.append(time)
            self.latitudes.append(float(observations[8]))
            self.longitudes.append(float(observations[7]))
            self.heights.append(float(observations[6]))
            count += 1

        self.times=np.asarray(times)
        
    def read_hypack_raw_file(self, fullpath):
        
        # This function will currently only function provided that there are GGA sentences in the records.
        # You may update the function to include other positioning messages as well, but this is 
        # outside the scope of the class
                        
        # Check the File's existence
        if os.path.exists(fullpath):
            self.data_path = fullpath
            print('Opening GNSS data file:' + fullpath)
        else:  # Raise a meaningful error
            raise RuntimeError('Unable to locate the input file' + fullpath)
            
        # Open, read and close the file
        hypack_file = open(fullpath)
        hypack_content = hypack_file.read()
        hypack_file.close    

        # Split the file in lines
        hypack_records = hypack_content.splitlines()
        
        # Go through the header lines to find the date of the survey (not contained in the GGA records)
        
        lines_parsed=0
        for hypack_record in hypack_records:
            
            # Check for the time and date
            
            if hypack_record[:3].lower() == "tnd":
                hypack_datetime=datetime.strptime(hypack_record[4:23], "%H:%M:%S %m/%d/%Y")
                
                print("HYPACK RAW Header start time and date: " + hypack_datetime.ctime())
                
            # Keep track of the lines parsed
            lines_parsed+=1

            # Stop going through the records if the record starts with the string eoh (End Of Header)
            if hypack_record[:3].lower() == "eoh":
                break         
        
        # We are at the end of the header - start looking for the first GGA record and compare its time 
        # to the TND record
        # This is so that we can set the correct date
        
        # Keep track of the number of GGA records found
        
        num_gga_recs=0
        
        for hypack_record in hypack_records[lines_parsed:]:

            if hypack_record[19:22] == "GGA":
                gga_data=hypack_record.split()[3]
                gga_data=gga_data.split(',')
                
                # Determine the time of day from both the header and the GGA string
     
                gga_timedelta=timedelta(hours=int(gga_data[1][0:2]), \
                                        minutes = int(gga_data[1][2:4]), \
                                        seconds = int(gga_data[1][4:6]))
            
                hypack_timedelta=timedelta(hours = hypack_datetime.hour, \
                                           minutes = hypack_datetime.minute, \
                                           seconds = hypack_datetime.second)
            
                # If the hours are not the same we need to believe one or the other, which one?
                # We know from experience that the GGA time from Trimble receivers is indeed utc, whereas we also know
                # that in older version of HYPACK the timestamp depended on the CPU time
                
                # Why do we care? We may end up with the timing for the wrong day if we are not careful! 
                # To know the correct date we need to take the date and time from the header and adjust it to UTC,
                # After that we know the time. If then the time difference between any succeeding GGA records 
                # is negative 
                # we know that we have just gone into a new day and need to update the datetime. 
                # Calculate the time difference in microseconds

                # Convert the HYPACK datetime to represent UTC time. We will make the ASSUMPTION that the header was 
                # created just before the first GGA record is received - this is _NOT_ fool proof (why?) 
                
                # On Piazza discuss a method to make this more robust by taking the time difference in minutes and 
                # seconds between the header time and the first GGA record time into account
                                
                hypack_datetime = hypack_datetime + \
                    timedelta( hours = round((gga_timedelta - hypack_timedelta).total_seconds() / 3600))
                
                print("HYPACK RAW Header start time and date in UTC: " + hypack_datetime.ctime())
                
                gga_datetime = hypack_datetime
                
                break
                
            # Keep track of where we are
                
            lines_parsed+=1
            
                
        # We are in the data file at the first GGA record - we know the date
  
        prev_gga_time = gga_timedelta
   
        
        for hypack_record in hypack_records[lines_parsed:]:
            if hypack_record[19:22] == "GGA":
                
                # Update the number of GGA records found
                
                num_gga_recs += 1
                
                # Get the GGA string and tokenize it
                
                gga_data=hypack_record.split()[3]
                gga_data=gga_data.split(',')
                
                # Determine the time of day from both the header and the GGA string
     
                gga_timedelta=timedelta(hours=int(gga_data[1][0:2]), \
                                        minutes = int(gga_data[1][2:4]), \
                                        seconds = int(gga_data[1][4:6]))
            
                # Check whether we rolled into a new day - we assume time is always increasing
                
                if gga_timedelta < prev_gga_time:
                    # We have reached the next day - update the date
                    
                    gga_datetime += timedelta(days=1)
                    
                    print( "Passed midnight, updating date to: "+str(gga_datetime.date()))
                    
                # Create the time
                
                time=gga_datetime.replace(hour=int(gga_data[1][0:2]), \
                                        minute = int(gga_data[1][2:4]), \
                                        second = int(gga_data[1][4:6]))

                # Add the time object to the list of times
                self.times.append( time)
                
                # Parse the latitude

                if gga_data[3].lower() == "n":
                    self.latitudes.append( float(gga_data[2][0:2])+float(gga_data[2][2:])/60.)
                else:
                    self.latitudes.append(-float(gga_data[2][0:2])-float(gga_data[2][2:])/60.)               

                # Parse the longitudes

                if gga_data[5].lower == "w":
                    self.longitudes.append( float(gga_data[4][0:3])+float(gga_data[4][3:])/60.)
                else:
                    self.longitudes.append(-float(gga_data[4][0:3])-float(gga_data[4][3:])/60.)               

                # Parse the GNSS Quality indicator
                
                self.qualities.append(int(gga_data[6]))

                # Parse the number of GNSS satellites used for the solution
                
                self.n_sats.append(int(gga_data[7]))

                # Parse the HDOP Quality indicator
                
                self.hdops.append(float(gga_data[8]))

                # Parse the orthometric height 
                
                self.heights.append(float(gga_data[9]))
                
                # Generate an error if the units of the orthometric height is not meters
                                   
                if gga_data[10].lower() != "m":
                    raise RuntimeError('Orthomeric height units are not meters!')                
                
                # Parse the geoid ellipsoid separation
                                   
                self.separations.append(float(gga_data[9]))
                                   
                if gga_data[12].lower() != "m":
                    raise RuntimeError('Orthomeric height units are not meters!') 
                    
                # If there is more data then parse it
                
                if gga_data[13] != "":
                    self.corr_ages.append(float(gga_data[13]))
                    self.corr_stations.append(float(gga_data[14]))
                    
                # For now, ignore the checksum (this would become a computer science assignment
                
                # Make sure to update the previous gga time
                    
                prev_gga_time=gga_timedelta

        # Set the reference ellipsoid to WGS84

        self.metadata["ellipsoid_name"] = "WGS84"
        self.metadata["geoid_name"] = "EGM08"
        self.metadata["height_relative_to"] = "geoid"

        # Let the user know how many GGA records there were in the file
        print("HYPACK RAW file contains: " + str(num_gga_recs) + " GGA records")

    def write(self):
        
        fullpath, _ = os.path.splitext(self.data_path)
        
        fullpath = fullpath + "_pos.txt"
       
        # Check the File's existence

        if os.path.exists(fullpath):
            # Let the user now we are overwriting an existing file
            print('Overwriting file: ' + fullpath)
        else:  # Let the user know we are writing to a file
            print('Writing to file: ' + fullpath)
            
        output_file = open(fullpath, mode="w")  # mode="w" to open the file in writing mode
        
        output_file.write('time latitude longitude\n')
        for i in range(0,len(self.times)):
            line_content=str(self.times[i]) + " %012.8f %013.8f\n" % (self.latitudes[i],self.longitudes[i])
            output_file.write(line_content)
    
    def write_hotlink(self, hotlink_path):
        
        fullpath, _ = os.path.splitext(self.data_path)
        
        fullpath = fullpath + "_pos.txt"
       
        # Check the File's existence

        if os.path.exists(fullpath):
            # Let the user now we are overwriting an existing file
            self.data_path = fullpath
            print('Overwriting file: ' + fullpath)
        else:  # Let the user know we are writing to a file
            print('Writing to file: ' + fullpath)
            
        output_file = open(fullpath, mode="w")  # mode="w" to open the file in writing mode
        
        # Write the header

        output_file.write('date time latitude longitude path\n')
        
        # Determine the duration associated to the positions in this object
        
        start_time = min(self.times)
        duration = max(self.times) - start_time
    
        # For each position write the date time longitude latitude path and fraction
        
        for i in range(0,len(self.times)):
            fraction = (self.times[i] - start_time) / duration
            line_content=str(self.times[i]) \
            + " %012.8f %013.8f %s?%f\n" % \
            (self.latitudes[i], self.longitudes[i], hotlink_path, fraction)
            output_file.write(line_content)
        
    def draw(self, projection='auto'):

        print('Drawing Positioning Data')

        # Dertermine the central latitude nand longitude
        central_lat = (max(self.latitudes)+min(self.latitudes))/2
        central_lon = (max(self.longitudes)+min(self.longitudes))/2

        # Create a 18x6 figure object
        fig = plt.figure(figsize=(18, 6))

        # Add a supertitle for the figure
        fig.suptitle("Positioning Data")

        # Create an Orthographic coordinate reference system (crs)
        crs_ortho = ccrs.Orthographic(central_lon, central_lat)

        # Plot the globe using the just defined crs in the first plot in a row of two subplots
        ax1 = fig.add_subplot(
            1, 2, 1, projection=crs_ortho)
        ax1.set_global()

        # Add Oceans, land and graticule
        ax1.add_feature(cartopy.feature.OCEAN)
        ax1.add_feature(cartopy.feature.LAND, edgecolor='black')
        ax1.gridlines()

        # Set the title
        ax1.set_title('Orthographic Map of Coverage Area')

        # Indicate the general area of the Positions by plotting a black marker with white edges
        # at the central latitude and longitude
        plt.plot(central_lon, central_lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                 markerfacecolor='black', markeredgecolor='white',
                 transform=crs_ortho)

        # Set the zone number and hemisphere for the UTM projection
        zone_number = int((np.floor((central_lon + 180) / 6) % 60) + 1)
        print(zone_number)
        if central_lat < 0:
            southern_hemisphere = True
        else:
            southern_hemisphere = False

        # Create an UTM reference system using the zone_number and hemisphere derived from the central coordinate
        crs_utm = ccrs.UTM(
            zone=zone_number, southern_hemisphere=southern_hemisphere)

        # Plot the zone using the UTM crs in the second plot of the two subplots
        ax2 = fig.add_subplot(
            1, 2, 2, projection=crs_utm)
        
        # Set display limits based on the extend of the geodetic coordinates in this object
        e_buffer=(np.max(self.longitudes)-np.min(self.longitudes))/10
        n_buffer=(np.max(self.latitudes)-np.min(self.latitudes))/5
        ax2.set_extent((np.min(self.longitudes)-e_buffer, np.max(self.longitudes)+e_buffer, np.min(
            self.latitudes)-n_buffer, np.max(self.latitudes)+n_buffer), crs=crs_utm)

        # Add Oceans, land and graticule
        ax2.add_feature(cartopy.feature.OCEAN, zorder=0)
        ax2.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        ax2.gridlines()
        
        # Set the title 
        ax2.set_title('UTM Map of Positioning Data')


        # Plot all the positions as black dots
        for lon, lat in zip(self.longitudes, self.latitudes):
            plt.plot(lon, lat, marker='.', markersize=2.0,
                     markerfacecolor='black', transform=crs_utm)

        # Display the figure
        plt.show()
        
        
    def carto_project(self, projection_name, z_reference):

        # start by creating a projection parameter string following PROJ protocol

        if not isinstance(projection_name, str):
            raise RuntimeError(
                'Position.carto_project(): argument `projection_name` must be of type str')

        if self.metadata["ellipsoid_name"] == None:
            raise RuntimeError(
                'Position.carto_project(): Requires ellipsoid metadata to be defined!')

        # Keep a list of projection that are implemented in this function
        implemented_projections = list()
        implemented_projections.append('utm')

        # Raise an error of a non implemented projection is asked for
        if projection_name.lower() not in implemented_projections:
            raise RuntimeError(
                'Position.project(): The projection `' + projection_name + '` is not yet implemented')
            
        # Universal will map to utm or ups depending on the central latitude
        if projection_name.lower() == 'universal':
            projection_name='utm'

        if projection_name.lower() == 'utm':
            proj_str = '+proj=utm'

            # Determine the central longitude and latitude in degrees
            central_lon = (max(self.longitudes)+min(self.longitudes))/2
            central_lat = (max(self.latitudes)+min(self.latitudes))/2

            # Determine the UTM Zone from the central longitude and latitude
            proj_str += ' +zone=' + \
                str(int((np.floor((central_lon + 180) / 6) % 60) + 1))
            if central_lat > 0:
                proj_str += ' +north'
            else:
                proj_str += ' +south'

            # The geodetic datum of the input coordinates
            proj_str += ' +ellps=' + self.metadata["ellipsoid_name"]

            # The datum for the output coordinates
            proj_str += ' +datum=' + self.metadata["ellipsoid_name"]

            # The units for the output coordinates
            proj_str += ' +units=m'

            # Prevent PROJ default behavior
            proj_str += ' +no_defs'

            # Create a pyproj object using proj_str
            proj_obj = proj.Proj(proj_str)

            # Perform the projection
            E, N = proj_obj(self.longitudes, self.latitudes)
            
            # Create a matrix of positions as 3D column vectors
            if z_reference.lower()=='ortho':
                self.proj_pos=np.asarray([np.asarray(E),np.asarray(N),np.asarray(self.heights)])
            else:
                raise RuntimeError('Position.carto_project(): currently only implemented for orthometric heights')
            
            # Add the string to the metadata
            self.metadata['proj_str'] = proj_str
            
            
            
            
def ParseNMEA0183_ZDA( dt_str):
    obs = dt_str.split(',')
    time = datetime( 
                int( obs[4]), int( obs[3]), int( obs[2]), 
                int( obs[1][0:2]), int( obs[1][2:4]), int( obs[1][4:6]),int(obs[1][7:])*10000)
    return time


def ParseNMEA0183_GGA( dt_str):
    
    # Get the GGA string and tokenize it
    gga_data = dt_str.split(',')

    # verify that we have a GGA string
    if not gga_data[0][-3:] == "GGA":
        raise RuntimeError(
                'ParseNMEA0183_GGA: argument `dt_str` must be a GGA message')

                
    # Determine the time of day from both the header and the GGA string
     
    gga_timedelta=timedelta(hours=int(gga_data[1][0:2]), \
                             minutes = int(gga_data[1][2:4]), \
                             seconds = int(gga_data[1][4:6]))
           
                
    # Parse the latitude
    if gga_data[3].lower() == "n":
        lat = float(gga_data[2][0:2])+float(gga_data[2][2:])/60.
    else:
        lat = -float(gga_data[2][0:2])-float(gga_data[2][2:])/60.             

    # Parse the longitude
    if gga_data[5].lower == "w":
        lon = float(gga_data[4][0:3])+float(gga_data[4][3:])/60.
    else:
        lon = -float(gga_data[4][0:3])-float(gga_data[4][3:])/60.               

    # Parse the GNSS Quality indicator
    q = int(gga_data[6])

    # Parse the number of GNSS satellites used for the solution
    n_sats = int(gga_data[7])

    # Parse the HDOP Quality indicator
    hdop = float(gga_data[8])

    # Parse the orthometric height 
    height = float(gga_data[9])
                
    # Generate an error if the units of the orthometric height is not meters
                                   
    if gga_data[10].lower() != "m":
        raise RuntimeError('Orthomeric height units are not meters!')                
                
    # Parse the geoid ellipsoid separation
    separation = float(gga_data[9])
                                   
    if gga_data[12].lower() != "m":
        raise RuntimeError('Orthomeric height units are not meters!') 
              
       
    # If there is more data then parse it
    corr_age = None
    corr_station = None
    if not gga_data[13] == "":
        corr_age = float(gga_data[13])
        corr_station = float(gga_data[14][0:-3])
                    
    # For now, ignore the checksum (this would become a computer science assignment

    return gga_timedelta, lat, lon, q, n_sats, hdop, height, separation, corr_age, corr_station