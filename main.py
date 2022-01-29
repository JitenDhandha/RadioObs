#Standard library components
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.constants as const
#Required astropy library components
from astropy import units as u
from astropy import coordinates as coords
from astropy import time as astrotime

######################################################################
#                   CHANGE THE VARIABLES BELOW                       #
#               Values given are only as an example                  #
######################################################################

#Enter the details of your telescope here
T_name = "42-ft telescope"
T_latitude = "53d 14m 10.5s"		#Latitude (degree,min,sec)
T_longitude = "2d 18m 25.7s"		#Longitude (degree,min,sec)
T_diameter = 12.8			#Diameter (metres)
T_frequency = 610			#Observing frequency (MHz)
T_bandwidth = 10			#Bandwidth (MHz)
T_aperture_efficiency = 0.55		#Aperture efficiency (dimensionless)
T_system_temperature = 130		#System temperature (Kelvin)

#Enter the details of the source you're trying to observe here
S_name = "B0833-45"					
S_RA = "08:35:20.6"			#Right Ascension (hour,min,sec)
S_Dec = "-45:10:34.8"			#Declination (degree,min,sec)
S_flux_density = 1100.00/1000 		#Flux density at telescope observing frequency (Jansky)


######################################################################
#                   	   Class: Telescope                          #
######################################################################

class Telescope:
	def __init__(self, t_name, t_location, t_diameter, t_frequency, t_bandwidth, t_aperture_efficiency, t_system_temperature):
		
		#Identifiers for the telescope
		self.name = t_name					#Type: str
		self.location = t_location				#Type: astropy.coordinates.EarthLocation, Units: Lat (degree,min,sec) and Lon (degree,min,sec)
		self.diameter = t_diameter				#Type: double, Units: m

		#Characteristics of the telescope
		self.frequency = t_frequency				#Type: double, Units: MHz
		self.bandwidth = t_bandwidth				#Type: double, Units: MHz
		self.aperture_efficiency = t_aperture_efficiency	#Type: double (between 0 and 1)
		self.system_temperature = t_system_temperature		#Type: double, Units: K

	#Effective area of the telescope (in m^2)
	def effective_area(self):

		effective_area = self.aperture_efficiency * const.pi * (self.diameter/2)**2
		return effective_area

	#Gain of the telescope (in K/Jy)
	def gain(self):

		gain = self.effective_area()/(2*const.k)
		gain *= 10**(-26)	#Conversion from K/(W/m^2/Hz) to K/Jy
		return gain

	#Function to print the details of the telescope
	def print_data(self):

		print("=========================================")
		print("TELESCOPE")
		print("Name: {0}".format(self.name))
		print("Lat / Lon: {0:.3f} / {1:.3f}".format(self.location.lat, self.location.lon))
		print("Diameter: {0} m".format(self.diameter))
		print("Observing frequency: {0} MHz".format(self.frequency))
		print("Bandwidth: {0} MHz".format(self.bandwidth))
		print("Aperture efficiency: {0} MHz".format(self.aperture_efficiency))
		print("System temperature: {0} K".format(self.system_temperature))
		print("=========================================")


######################################################################
#                   	     Class: Source                           #
######################################################################

class Source:
	def __init__(self, s_name, s_location, s_flux_density):

		self.name = s_name 					#Type: str
		self.location = s_location 				#Type: astropy.coordinates.SkyCoord, Units: RA (hour,min,sec) and Dec (degree,min,sec)
		self.flux_density = s_flux_density			#Type: double, Units: Jy

	#Function to print the details of the telescope
	def print_data(self):

		print("=========================================")
		print("SOURCE")
		print("Name: {0}".format(self.name))
		print("RA / Dec: {0:.3f} / {1:.3f}".format(self.location.ra, self.location.dec))
		print("Flux density (at telescope observing frequency): {0} Jy".format(self.flux_density))
		print("=========================================")

	#Function to plot the altitude of the source in the sky at a given location and day
	def plot_altitude(self, telescope, date):

		#Filling up array for different observation times
		dt = 30*60 	#in seconds
		dt_astropy = astrotime.TimeDelta(dt, format='sec')
		observation_times = astrotime.Time(date + ' ' + '00:00:00', location=telescope.location) + dt_astropy * np.arange(0,86400//dt)

		#Calculating the altitude of the source at the different observation times
		altitudes = []
		for observation_time in observation_times:
			aa = coords.AltAz(location=telescope.location, obstime=observation_time)
			altitudes.append(self.location.transform_to(aa).alt.deg)

		#Plotting the altitudes as a function of time of day
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot_date(observation_times.plot_date, altitudes, fmt='ko', ls='-', label=r'$\mathrm{Source}$')
		ax.axhline(0,c='red',linestyle='-',label=r'$\mathrm{Horizon}$')
		ax.axhline(10,c='red',linestyle='--',label=r'$10^{\circ}~\mathrm{degree}$')
		ax.set_xlabel(r"$\mathrm{Time~of~the~day~[UTC,~hours]}$",fontsize='large')
		ax.set_ylabel(r"$\mathrm{Altitude~in~the~sky~[degrees]}$",fontsize='large')
		ax.set_title(r"$\mathrm{{Altitude~of~source~{0}~in~the~sky~on~{1}}}$".format(self.name,date),fontsize='large')
		fig.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))	#Time axis formatting
		ax.legend()
		fig.show()

	#Function to print the integration time for the telescope observation
	def print_integration_time(self, telescope, SNR):

		#Telescope bandwidth converted from MHz to Hz.
		#The 2 in the denominator below denotes the number of polarizations.
		integration_time = SNR**2 * (telescope.system_temperature**2) / ((telescope.gain() * (self.flux_density))**2 * 2 * telescope.bandwidth*10**(6))
		print("The integration time for observing at SNR of {0:.1f} is {1:.3f} seconds.".format(SNR, integration_time))


######################################################################
#                   	       Main code                             #
######################################################################

def main():

	#Creating object of type Telescope and printing details
	T_location = coords.EarthLocation(lat=T_latitude, lon=T_longitude)
	Tel = Telescope(T_name, T_location, T_diameter, T_frequency, T_bandwidth, T_aperture_efficiency, T_system_temperature)
	Tel.print_data()

	#Creating object of type Source and printing details
	S_location = coords.SkyCoord(ra=S_RA, dec=S_Dec, unit=(u.hourangle, u.deg))
	Src = Source(S_name, S_location, S_flux_density)
	Src.print_data()

	#Taking inputs and printing outputs
	doplot = input("Would you like to plot the altitude of the source in the sky? ")
	if(doplot.lower() in ["y","yes"]):
		date = input("Enter date of observation in YYYY-MM-DD format: ")
		try:
			astrotime.Time(date + ' ' + '00:00:00')
			Src.plot_altitude(Tel, date)
		except ValueError:
			print("Incorrect format for date!")
			return
	
	dotime = input("Would you like to know your integration time? ")
	if(dotime.lower() in ["y","yes"]):
		SNR = input("Enter signal-to-noise ratio you would like to achieve: ")
		try:
			SNR = float(SNR)
			Src.print_integration_time(Tel, SNR)
		except ValueError:
			print("Not a valid signal-to-noise ratio!")
			return

if __name__ == "__main__":
	main()
