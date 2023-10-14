import os.path
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

from bisect import bisect_left
from scipy.interpolate import make_interp_spline, BSpline

	
# parse flight sequence data
def read_flight_data(fileName, sh=1):
	ldata = np.genfromtxt(fileName, delimiter='|', skip_header=sh, autostrip=True, dtype='unicode')
	# split the arrays
	event = ldata[:,0]
	time = ldata[:,1]
	altitude = ldata[:,2]
	velocity = ldata[:,3]
	# strip str
	event = [i.strip() for i in event]
	# convert str to float	
	time = [float(i) for i in time]
	altitude = [float(i) for i in altitude]
	velocity = [float(i) for i in velocity]	
	return event, time, altitude, velocity

# parse velocity data
def read_velocity_data(fileName, sh=1):
	ldata = np.genfromtxt(fileName, delimiter='|', skip_header=sh, autostrip=True, dtype='unicode')
	# split the arrays
	time = ldata[:,0]
	velocity = ldata[:,1]
	# convert str to float
	time = [float(i) for i in time]
	velocity = [float(i) for i in velocity]	
	return time, velocity

# parse altitude data
def read_altitude_data(fileName, sh=1):
	ldata = np.genfromtxt(fileName, delimiter='|', skip_header=sh, autostrip=True, dtype='unicode')
	# split the arrays
	time = ldata[:,0]
	altitude = ldata[:,1]
	# convert str to float
	time = [float(i) for i in time]
	altitude = [float(i) for i in altitude]	
	return time, altitude

# parse atmospheric data
def read_atmospheric_data(fileName='./atm_data.dat', sh=2):
	atw = np.genfromtxt(fileName, delimiter='|', skip_header=sh)
	return atw
	
# simple linear interpolation
def interpolate(x_list, y_list, x):  
	if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
		raise ValueError("x_list must be in strictly ascending order!")
	intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
	slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals] 

	if x <= x_list[0]:
		return y_list[0]
	elif x >= x_list[-1]:
		return y_list[-1]
	else:
		i = bisect_left(x_list, x) - 1
		return y_list[i] + slopes[i] * (x - x_list[i])

# returns interpolated atmos data for each alt (altitude value)	
def get_atmospheric_data(ndarray, alt):
	# split the arrays
	altitude = ndarray[:,0]
	temperature = ndarray[:,1]
	density = ndarray[:,2]
	pressure = ndarray[:,3]
	molwt = ndarray[:,4]
	# interpolate values
	t = interpolate(altitude, temperature, alt)
	d = interpolate(altitude, density, alt)
	p = interpolate(altitude, pressure, alt)
	m = interpolate(altitude, molwt, alt)	
	return t, d, p, m

# find the nearest index in the list for the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

# get atmospheric data for each altitude. Units returned (kg/m3)
def get_density_data(file_atm, altitude):	
	temperature = []; pressure = []; density = [];
	# file check
	if os.path.exists(file_atm):
		# read Altitude(km), Temp(K), Density(kg/m3), Pressure(Pa), Mol. Wt.(kg/kmol)
		ndr = read_atmospheric_data(file_atm, 2)	# skip first 2 headers
		# get atmos data for each altitude 
		for index, elem in enumerate(altitude):
			t, d, p, m = get_atmospheric_data(ndr, elem )
			density.append(d)
	return density

# calculate avg acceleration from velocity and time intervel
def get_acceleration(velocity, time):
	acceleration = []
	pre_time = 0; pre_velocity = 0;
	for i, (t, v) in enumerate( zip(time, velocity)): 
		dv = v - pre_velocity
		dt = t - pre_time			
		if dt != 0:  a = dv/dt
		else: a = 0
		pre_time = t; pre_velocity = v;		
		acceleration.append(a)	
	return acceleration

# calculate avg jerk from acceleration and time intervel
def get_jerk(acceleration, time):
	jerk = []
	pre_time = 0; pre_acceleration = 0;
	for i, (t, a) in enumerate( zip(time, acceleration)): 
		da = a - pre_acceleration
		dt = t - pre_time			
		if dt != 0:  j = da/dt
		else: j = 0
		pre_time = t; pre_acceleration = a;		
		jerk.append(j)	
	return jerk

# set font sizes
def set_plot_fonts():
	# font size
	S_SIZE = 9
	plt.rc('font', size=S_SIZE)         # controls default text sizes
	plt.rc('axes', titlesize=S_SIZE)    # fontsize of the axes title
	plt.rc('axes', labelsize=S_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=S_SIZE)   # fontsize of the tick labels
	plt.rc('ytick', labelsize=S_SIZE)   # fontsize of the tick labels	 
	return

# plot altitude and relative velocity vs time
def plot_merged_display(plot_label, time, altitude, velocity, events_data):
	fig, ax1 = plt.subplots( figsize=(12,5))
	fig.suptitle(plot_label, fontsize=10)
	fig.subplots_adjust(top=0.8)	
	# fix font sizes
	set_plot_fonts()
	S_SIZE = 8
	plot_text = 'Merged_Display'
	ax1.set_title(plot_text, fontsize=S_SIZE)
	# set the x label
	ax1.set_xlabel('Time (sec)')
	ax1.set_ylabel('Altitude (m)', color="red")
	ax1.plot(time, altitude, color='red')
	ax1.tick_params(axis='y', colors='red')
	# instantiate a second axes that shares the same x-axis
	ax2 = ax1.twinx()  
	# we already handled the x-label with ax1
	ax2.set_ylabel('Relative Velocity (km/s)', color="blue")
	ax2.plot(time, velocity, color='blue')
	ax2.tick_params(axis='y', colors='blue')
	# mark events. events_data = (event, time, altitude, velocity)
	for i, (event, px)  in enumerate(zip(events_data[0], events_data[1])):
		time_index, val = find_nearest(time, px)
		py1 = altitude[time_index]
		py2 = velocity[time_index]
		ax1.scatter( px, py1, s=30, color='g')
		ax2.scatter( px, py2, s=30, color='g')
		ax1.text( px, py1, event, fontsize='x-small')
		ax2.text( px, py2, event, fontsize='x-small')
	# color
	plot_bg_color = (1, 1, 0.5)	
	# draw ticks like a graph sheet
	ax1.set_facecolor(plot_bg_color)	
	# Turn on the minor TICKS, which are required for the minor GRID
	ax1.minorticks_on()
	# Customize the major grid
	ax1.grid(which='major', linestyle='-', linewidth='0.5', color='red')
	# Customize the minor grid
	ax1.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
	# Turn off the display of all ticks.
	ax1.tick_params(which='both', # Options for both major and minor ticks
					top='off', # turn off top ticks
					left='off', # turn off left ticks
					right='off',  # turn off right ticks
					bottom='off') # turn off bottom ticks	
	
	# show
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	s_name = plot_label.lower() + '_' + plot_text.lower() + '.png'
	plt.savefig(s_name)
	plt.show()
	return

# plot jerk vs time	
def plot_jerk(plot_label, time, jerk, events_data):
	# Create a figure with 2 rows and 2 cols of subplots
	fig, ax = plt.subplots( figsize=(12,5) )
	fig.suptitle(plot_label, fontsize=10)
	fig.subplots_adjust(top=0.8)	
	plot_text = 'Jerk_vs_Time'
	alt_label = "Jerk (in G/s)"
	time_label = "Time (sec)"
	# fix font sizes
	set_plot_fonts()
	S_SIZE = 8
	ax.set_title(plot_text, fontsize=S_SIZE)
	ax.tick_params(axis='both', which='major', labelsize=S_SIZE) 
	ax.tick_params(axis='both', which='minor', labelsize=S_SIZE)	
	# color
	plot_bg_color = (1, 1, 0.5)
	# plot
	ax.plot(time, jerk)
	# mark events. events_data = (event, time, altitude, velocity)
	for i, (event, px)  in enumerate(zip(events_data[0], events_data[1])):
		time_index, val = find_nearest(time, px)
		py = jerk[time_index]
		ax.scatter( px, py, s=30, color='g')
		ax.text( px, py, event, fontsize='x-small')
	# set axis labels 	
	# ax.xaxis.set_label_coords(0.5, -0.06)
	ax.set_xlabel(time_label, fontsize=S_SIZE)
	# ax.yaxis.set_label_coords(-0.1,0.5)
	ax.set_ylabel(alt_label, fontsize=S_SIZE)
	# draw ticks like a graph sheet
	ax.set_facecolor(plot_bg_color)
	# turn on the minor TICKS
	ax.minorticks_on()
	# customize the major grid
	ax.grid(which='major', linestyle='-', linewidth='0.5', color='red')
	# customize the minor grid
	ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')		

	# show
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	s_name = plot_label.lower() + '_' + plot_text.lower() + '.png'
	plt.savefig(s_name)
	plt.show()
	return


# plot acceleration vs time	
def plot_acceleration(plot_label, time, acceleration, events_data):
	# Create a figure with 2 rows and 2 cols of subplots
	fig, ax = plt.subplots( figsize=(12,5) )
	fig.suptitle(plot_label, fontsize=10)
	fig.subplots_adjust(top=0.8)	
	plot_text = 'Acel_vs_Time'
	alt_label = "Accleration (in G's)"
	time_label = "Time (sec)"
	# fix font sizes
	set_plot_fonts()
	S_SIZE = 8
	ax.set_title(plot_text, fontsize=S_SIZE)
	ax.tick_params(axis='both', which='major', labelsize=S_SIZE) 
	ax.tick_params(axis='both', which='minor', labelsize=S_SIZE)	
	# color
	plot_bg_color = (1, 1, 0.5)
	# plot
	ax.plot(time, acceleration)
	# mark events. events_data = (event, time, altitude, velocity)
	for i, (event, px)  in enumerate(zip(events_data[0], events_data[1])):
		time_index, val = find_nearest(time, px)
		py = acceleration[time_index]
		ax.scatter( px, py, s=30, color='g')
		ax.text( px, py, event, fontsize='x-small')
	# set axis labels 	
	# ax.xaxis.set_label_coords(0.5, -0.06)
	ax.set_xlabel(time_label, fontsize=S_SIZE)
	# ax.yaxis.set_label_coords(-0.1,0.5)
	ax.set_ylabel(alt_label, fontsize=S_SIZE)
	# draw ticks like a graph sheet
	ax.set_facecolor(plot_bg_color)
	# turn on the minor TICKS
	ax.minorticks_on()
	# customize the major grid
	ax.grid(which='major', linestyle='-', linewidth='0.5', color='red')
	# customize the minor grid
	ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')		

	# show
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	s_name = plot_label.lower() + '_' + plot_text.lower() + '.png'
	plt.savefig(s_name)
	plt.show()
	return

# plot dynamic pressure vs alt
def plot_dynamic_pressure(plot_label, altitude, dynamic_pressure, events_data):
	# Create a figure with 2 rows and 2 cols of subplots
	fig, ax = plt.subplots( figsize=(12,5) )
	fig.suptitle(plot_label, fontsize=10)
	fig.subplots_adjust(top=0.8)	
	plot_text = 'Alt_vs_Dyn_Pressure'
	x_label = "Altitude (km)"
	y_label = "Dynamic Pressure (kPa)"	
	# setup fonts
	set_plot_fonts()
	S_SIZE = 8
	ax.set_title(plot_text, fontsize=S_SIZE)
	ax.tick_params(axis='both', which='major', labelsize=S_SIZE) 
	ax.tick_params(axis='both', which='minor', labelsize=S_SIZE)	
	# color
	plot_bg_color = (1, 1, 0.5)
	# plot dynamic pressure vs altitude
	ax.plot(altitude, dynamic_pressure)
	# mark events. events_data = (event, time, altitude, velocity)
	for i, (event, px)  in enumerate(zip(events_data[0], events_data[2])):
		py = 0
		ax.scatter( px, py, s=30, color='g')
		ax.text( px, py+1, event, rotation=90, va='bottom', fontsize='x-small')
	# find max-q
	ymax = max(dynamic_pressure)
	xpos = dynamic_pressure.index(ymax)
	xmax = altitude[xpos]
	max_label = 'Max-Q = ' + str(round(ymax,2)) + ' kPa at altitude '+ str(round(xmax,2)) + ' km '
	# indicate max-q
	ax.annotate(max_label, xy=(xmax, ymax), xytext=(xmax, ymax-10),
				arrowprops=dict(facecolor='black', shrink=0.1, width=0.5) )
	# set axis labels
	# ax.xaxis.set_label_coords(0.5, -0.06)
	ax.set_xlabel(x_label, fontsize=S_SIZE)
	# ax.yaxis.set_label_coords(-0.1,0.5)
	ax.set_ylabel(y_label, fontsize=S_SIZE)
	# draw ticks like a graph sheet
	ax.set_facecolor(plot_bg_color)
	# turn on the minor TICKS
	ax.minorticks_on()
	# customize the major grid
	ax.grid(which='major', linestyle='-', linewidth='0.5', color='red')
	# customize the minor grid
	ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')	
	# show
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	s_name = plot_label.lower() + '_' + plot_text.lower() + '.png'
	plt.savefig(s_name)
	plt.show()
	return

# butterworth filter
def butterworth(rough_data):
	# Filter order
	N  = 7 
	# Cutoff frequency
	Wn = 0.07 
	B, A = signal.butter(N, Wn, output='ba')
	smooth_data = signal.filtfilt(B,A, rough_data)
	return smooth_data
    	    
# generate data for plotting acceleration vs time
def get_acceleration_plot_data (file_evt, file_vel, low_noise_filter_flag = False ):
	# file check
	if os.path.exists(file_evt):
		# read event, time, altitude, velocity
		event, time, altitude, velocity = read_flight_data(file_evt, 1)
		events_data = (event, time, altitude, velocity)		
		# relative velocity
		velocity = [x - velocity[0] for x in velocity]		
	# read velocity data (sec, m/sec)
	f_time, f_velocity = read_velocity_data(file_vel, 1)
	# f_time = [round(t) for t in f_time]
	# get avg acceleration
	f_acceleration = get_acceleration(f_velocity, f_time)	
	# round it for plotting
	f_time = [round(t) for t in f_time]	
	# acceleration (m/s2) to (in G's) for plotting
	f_acceleration = [a/9.8 for a in f_acceleration]
	# output the data (filtered vs non-filtered) based on the flag
	if low_noise_filter_flag: 
		# butterworth filter		
		f_acceleration_filtered = butterworth( np.array(f_acceleration))			
		accl_data = (f_time, f_acceleration_filtered.tolist())
	else: 
		accl_data = (f_time, f_acceleration)	
	# accleration data (f_time, f_acceleration data's)
	return accl_data, events_data

# generate data for plotting acceleration vs time
def get_jerk_plot_data (accl_data, events_data, low_noise_filter_flag = False ):
    (f_time, f_acceleration) = accl_data
    (event, time, altitude, velocity) = events_data
    # calculate avg jerk from acceleration and time intervel
    f_jerk =  get_jerk(f_acceleration, f_time)
    # acceleration (m/s2) to (in G's) for plotting
    #f_acceleration = [a/9.8 for a in f_acceleration]
    # output the data (filtered vs non-filtered) based on the flag
    if low_noise_filter_flag: 
        # butterworth filter		
        f_jerk_filtered = butterworth( np.array(f_jerk))			
        jerk_data = (f_time, f_jerk_filtered.tolist())
    else: 
        jerk_data = (f_time, f_jerk)	
    # accleration data (f_time, f_acceleration data's)
    return jerk_data, events_data


# generate data for plotting dynamic pressure vs alt	
def get_dynamic_pressure_plot_data ( file_evt, file_alt, file_vel, filename_atm, low_noise_filter_flag = False ):
	combined_time = []; combined_velocity = []; combined_altitude = [];
	# file check
	if os.path.exists(file_evt):
		# read event, time, altitude, velocity
		event, time, altitude, velocity = read_flight_data(file_evt, 1)
		events_data = (event, time, altitude, velocity)		
		# relative velocity
		velocity = [x - velocity[0] for x in velocity]		
	# read velocity data (sec, m/sec)
	f_time, f_velocity = read_velocity_data(file_vel, 1)
	# read altitude data (sec, km)
	h_time, h_altitude = read_altitude_data(file_alt, 1)	
	# combine data
	f_time = [round(t) for t in f_time]	
	h_time = [round(t) for t in h_time]
	for index, (ht, ha) in enumerate(zip(h_time, h_altitude)):
		if ht in f_time:
			f_index = f_time.index(ht)
			combined_time.append(ht)
			combined_altitude.append(ha)
			combined_velocity.append(f_velocity[f_index]) 
	# get atmospheric data
	combined_density = get_density_data(filename_atm, combined_altitude)
	# variables
	f_dynamic = [];
	# calculate dynamic pressure and drag force
	for a, v, d in zip(combined_altitude, combined_velocity, combined_density):
		q  = 1/2 * d * v**2;	# dynamic pressure
		f_dynamic.append(q)

	# velocity (m/s) to (Km/s) for plotting
	combined_velocity = [v/1000 for v in combined_velocity]
	# Pascal to KPa for plotting
	f_dynamic = [p/1000 for p in f_dynamic]	

	# tavd_data (time(s), altidude(Km), velocity(Km/s), dynamic pressure(KPa) data's)
	tavd_data = (combined_time, combined_altitude, combined_velocity, f_dynamic)	
	# output the data (filtered vs non-filtered) based on the flag
	if low_noise_filter_flag: 
		# butterworth filter		
		f_dynamic_filtered = butterworth( np.array(f_dynamic))
		# tavd_data (time(s), altidude(Km), velocity(Km/s), dynamic pressure(KPa) data's)
		tavd_data = (combined_time, combined_altitude, combined_velocity, f_dynamic_filtered.tolist())		
	else: 
		# tavd_data (time(s), altidude(Km), velocity(Km/s), dynamic pressure(KPa) data's)
		tavd_data = (combined_time, combined_altitude, combined_velocity, f_dynamic)	
	
	return tavd_data, events_data

# plot vehicle profiles
def plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag ):
	# get acceleration data
	accl_data, events_data = get_acceleration_plot_data (filename_evt, filename_vel, apply_low_noise_filter_flag )
	# plot acceleration
	plot_acceleration(vehicle_mission_name, accl_data[0], accl_data[1], events_data)
	# get jerk data (time, jerk)
	jerk_data, events_data = get_jerk_plot_data (accl_data, events_data, apply_low_noise_filter_flag )
	# plot jerk
	plot_jerk(vehicle_mission_name, jerk_data[0], jerk_data[1], events_data)	
	# get time, altitude, velocity, dynamic pressure data
	tavd_data, events_data = get_dynamic_pressure_plot_data ( filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )
	# plot dynamic pressure
	plot_dynamic_pressure(vehicle_mission_name, tavd_data[1], tavd_data[3], events_data)
	# plot merged display
	plot_merged_display(vehicle_mission_name, tavd_data[0], tavd_data[1], tavd_data[2], events_data)

	return
	
# main function
if __name__ == "__main__":
	# atmospheric data
	filename_atm = './data/atm_data.dat'
	
	# plot data with noise filtering 
	apply_low_noise_filter_flag = True

	# GSLV-MK3-D2 Data
	vehicle_mission_name = 'GSLV-MK3-D2-GSAT-29'
	filename_evt = './data/gslv-mk3-d2-gsat29-flight-events.dat'
	filename_alt = './data/gslv-mk3-d2-gsat29-48-alt.dat'
	filename_vel = './data/gslv-mk3-d2-gsat29-48-vel.dat'
	# plot gslv-mk3 profiles	
	plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )

<<<<<<< HEAD
	# # GSLV-MK2-F08 Data
	# vehicle_mission_name = 'GSLV-MK2-F08-GSAT-6A'
	# filename_evt = './data/gslv-f08-gsat6a-flight-events.dat'
	# filename_alt = './data/gslv-f08-gsat6a-42-alt.dat'
	# filename_vel = './data/gslv-f08-gsat6a-42-vel.dat'	
	# # plot gslv-mk2 profiles	
	# plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )

	# # PSLV-C26 Data
	# vehicle_mission_name = 'PSLV-C26-IRNSS-1C'
	# filename_evt = './data/pslv-c26-irnss-1c-flight-events.dat'
	# filename_alt = './data/pslv-c26-irnss-1c-11-alt.dat'
	# filename_vel = './data/pslv-c26-irnss-1c-11-vel.dat'
	# # plot pslv profiles
	# plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )
	
	vehicle_mission_name = 'GSLV-MK3-M4-C3'
	filename_evt = './data/gslv-mk3-m4-c3-flight-eventsMod.dat'
	filename_alt = './data/gslv-mk3-m4-c3-alt.dat'
	filename_vel = './data/gslv-mk3-m4-c3-vel.dat'	
	# plot gslv-mk2 profiles	
=======
	# GSLV-MK2-F08 Data
	vehicle_mission_name = 'GSLV-MK2-F08-GSAT-6A'
	filename_evt = './data/gslv-f08-gsat6a-flight-events.dat'
	filename_alt = './data/gslv-f08-gsat6a-42-alt.dat'
	filename_vel = './data/gslv-f08-gsat6a-42-vel.dat'	
	# plot gslv-mk2 profiles	
	plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )

	# PSLV-C26 Data
	vehicle_mission_name = 'PSLV-C26-IRNSS-1C'
	filename_evt = './data/pslv-c26-irnss-1c-flight-events.dat'
	filename_alt = './data/pslv-c26-irnss-1c-11-alt.dat'
	filename_vel = './data/pslv-c26-irnss-1c-11-vel.dat'
	# plot pslv profiles
>>>>>>> f9c7007050b3c6188b5e16a9c3e782dea9727284
	plot_profiles( vehicle_mission_name, filename_evt, filename_alt, filename_vel, filename_atm, apply_low_noise_filter_flag )
