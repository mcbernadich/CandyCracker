import numpy as np
from scipy.optimize import newton
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import argparse
import sys

def deg_min_sec_to_deg(value):
	value=value.split(":")
	if len(value)==3:
		degrees=float(value[0])+float(value[1])/60+float(value[2])/3600
	elif len(value)==2:
		degrees=float(value[0])+float(value[1])/60
	elif len(value)==1:
		degrees=float(value[0])
	else:
		sys.exit("Declination must come in degrees, dd:mm or dd:mm:ss")
	return degrees

def h_min_sec_to_deg(value):
	value=value.split(":")
	if len(value)==3:
		degrees=180*(float(value[0])+float(value[1])/60+float(value[2])/3600)/12
	elif len(value)==1:
		degrees=float(value[0])
	else:
		sys.exit("Righ ascension must come in degrees or hh:mm:ss")
	return degrees

def phase_to_eccentric(phase,e,om):
	true=phase-om
	eccentric=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(true/2))
	return eccentric

def phase_to_eccentric_with_omdot(phase,e,om,omdot,t,t0):
	true=phase-om-omdot*(t-t0)
	eccentric=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(true/2))
	return eccentric

def passage_function_with_omdot(t,p,e,om,omdot,t0,k):
	eccentric_anomaly=phase_to_eccentric_with_omdot(phase,e,om,omdot,t,t0)
	return t-(t0+(p/2/np.pi)*(2*k*np.pi+eccentric_anomaly-e*np.sin(eccentric_anomaly)))

def get_times(phase,p,e,om,t0,k):
	eccentric_anomaly=phase_to_eccentric(phase,e,om)
	time=t0+(p/2/np.pi)*(2*k*np.pi+eccentric_anomaly-e*np.sin(eccentric_anomaly))
#	print(time)
	return time

def get_times_with_omdot(phase,p,e,om,omdot,t0,k):
	# First, get the trial.
	eccentric_anomaly=phase_to_eccentric(phase,e,om)
	trial_time=t0+(p/2/np.pi)*(2*k*np.pi+eccentric_anomaly-e*np.sin(eccentric_anomaly))
	# Then, run the rest numerically.
	time=newton(passage_function_with_omdot,trial_time,args=(p,e,om,omdot,t0,k))
#	print(trial_time,time)
	return time

def get_times_from_mean(mean_anomaly,p,t0,k):
	periastron_passage=t0+(k+mean_anomaly)*p
	return periastron_passage

def read_ephemeris(file):
	ephemeris=open(file,"r")
	ecc=0
	omega=0
	for line in ephemeris:
		line=line.strip().split()
		if line[0]=="PB":
			p_orb=float(line[1])
		elif line[0]=="ECC":
			ecc=float(line[1])
		elif line[0]=="OM":
			omega=float(line[1])
		elif line[0]=="T0":
			periastron=float(line[1])
		elif line[0]=="OMDOT":
			omdot=float(line[1])
		elif line[0]=="RAJ":
			ra=str(line[1])
		elif line[0]=="DECJ":
			dec=str(line[1])
	if f0 and p_orb and periastron and ecc and omega and omdot:
		print("Pulsar parameters loaded from {}:".format(file))
		print("- Righ ascension: "+ra)
		print("- Declination: "+dec)
		print("- Orbital period: {} days".format(p_orb))
		print("- Eccentricity: {}".format(ecc))
		print("- Periastron angle: {} degrees".format(omega*(np.pi/180)))
		print("- Periastron angle: {} º/yr".format(omdot*(np.pi/180)*(1/365)))
		print("- Periastron passage: {} (MJD)".format(periastron))
		return ra,dec,f0,p_orb,x,ecc,omega,periastron,omdot
	if f0 and p_orb and periastron and ecc and omega:
		print("Pulsar parameters loaded from {}:".format(file))
		print("- Righ ascension: "+ra)
		print("- Declination: "+dec)
		print("- Orbital period: {} days".format(p_orb))
		print("- Eccentricity: {}".format(ecc))
		print("- Periastron angle: {} degrees".format(omega*(np.pi/180)))
		print("- Periastron passage: {} (MJD)".format(periastron))
		return ra,dec,f0,p_orb,x,ecc,omega,periastron,"none"
	else:
		sys.exit("The ephemeris file doesn't include all the necessary parameters.")

parser=argparse.ArgumentParser(description="Take in an orbital model an compute OBSERVABLE superior conjunction and periastron passage times.")
parser.add_argument("--ephemeris",help="Fitorbit or tempo or tempo2-format ephemeris.")
parser.add_argument("--phase",type=float,help="Orbital phase for which times are computed, from 0 to 1.")
parser.add_argument("--phase_def",type=str,help="Phase can be either defined as angle from ascending node (true, from true anomaly) or mean anomaly (mean).", choices=["true","mean"])
parser.add_argument("--ra",type=str,help="Righ ascension in degrees or hh:mm:ss")
parser.add_argument("--dec",type=str,help="Declination in degrees, dd:mm or dd:mm:ss")
#parser.add_argument("--telescope",type=str,help="Location of observer", choices=["Effelsberg","MeerKAT","Parkes"])
parser.add_argument("-p","--period",type=float,help="Orbital period in days")
parser.add_argument("-e","--eccentricity",type=float,help="Excentricity. If not given, assumed 0.")
parser.add_argument("-o","--omega",type=float,help="Periastron angle in degrees. If not given, assumed 0.")
parser.add_argument("-t","--t_periastron",type=float,help="Periastron passage time/epoch of periastron, in mjd.")
parser.add_argument("-r","--range",help="Interval of integer orbits from the first periastron passage to output in the result, as in 'min:max'. Default: '0:10'.")
parser.add_argument("--tolerance",help="Minutes before and after the event that needs to be included in the observation as in 'before:after'. Default: '30:30'. Negative numbers accepted.")
parser.add_argument("--delay",type=float,help="Days to shift the center of the observation forward (positive) or backward (negative) in time.")
parser.add_argument("--omdot",type=float,help="Periastron advance in º/yr. Optional.")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

print(" ")

if args.phase or args.phase==0.0:
	phase=args.phase-int(args.phase)
	if args.verbose==True:
		print("Looking for observing times at phase "+str(phase))
		print("")
else:
	sys.exit("Please specify the desired orbital phase with --phase")

if args.phase_def:
	phase_def=args.phase_def
	if args.verbose==True:
		print("Orbital phase is assumed to be based on "+str(phase_def)+" anomaly.")
		if phase_def=="true":
			print("The required phase is counted from the ascending node.")
		if phase_def=="true":
			print("The required phase is counted from periastron.")
		print(" ")
else:
	phase_def="true"
	print("No phase measure definition. Assuming that true anomaly angle from ascending node is required.")
	print(" ")

if args.delay or args.delay==0.0:
	delay=args.delay*u.day
	if args.verbose==True:
		print("Introducing a delay of "+str(args.delay)+" days.")
		print("")
else:
	delay=0.0
	if args.verbose==True:
		print("No delay has been intorduced.")
		print("")

if args.tolerance:
	prior=float(args.tolerance.split(":")[0])/60
	posterior=float(args.tolerance.split(":")[1])/60
	if args.verbose==True:
		print("Planning observations for {} and {} minutes after and before the event.".format(prior*60,posterior*60))
		print("")
else:
	prior=0.5
	posterior=0.5
	if args.verbose==True:
		print("Planning observations for {} and {} minutes after and before the event.".format(prior*60,posterior*60))
		print("")

if args.ephemeris:
	(ra,dec,f0,p_orb,x,ecc,omega,periastron,omdot)=read_ephemeris(args.ephemeris)
	ra=h_min_sec_to_deg(ra)
	dec=deg_min_sec_to_deg(dec)
	if type(omdot)=="str":
		fit_omdot=False
	else:
		fit_omdot=True
else:

	if args.ra:
		ra=h_min_sec_to_deg(args.ra)
		if args.verbose==True:
			print("Righ ascension: "+str(ra))
			print("")
	else:
		sys.exit("Please specify right ascension with --ra")

	if args.dec:
		dec=deg_min_sec_to_deg(args.dec)
		if args.verbose==True:
			print("Declination: "+str(dec))
			print("")
	else:
		sys.exit("Please specify declination with --dec")

	if args.range:
		start=int(args.range.split(":")[0])
		end=int(args.range.split(":")[1])
		if start>end:
			end2=end
			end=start
			start=end2
		if args.verbose==True:
			print("Giving times from orbit {} {}.".format(start,end))
			print("")
	else:
		start=0
		end=10
		if args.verbose==True:
			print("Giving times from orbit {} {}.".format(start,end))
			print("")

	if args.period:
		p=args.period
		if args.verbose==True:
			print("Orbital period: {} days".format(p))
			print(" ")
	else:
		sys.exit("Please specify orbital period in days with -p")

	if args.t_periastron or args.t_periastron==0.0:
		t0=args.t_periastron
		if args.verbose==True:
			print("Periastron passage time: {} mjd".format(t0))
			print(" ")
	else:
		sys.exit("Please specify periastron passage time in mjd with --t0")

	if args.eccentricity:
		e=args.eccentricity
		if args.verbose==True:
			print("Eccentricity: {} ls".format(e))
			print(" ")
	else:
		e=0
		if args.verbose==True:
			print("Eccentricity assumed to be 0")
			print(" ")

	if args.omega:
		om=args.omega*(np.pi/180)
		if args.verbose==True:
			print("Periastron angle: {} degrees".format(om))
			print(" ")
	else:
		om=0
		if args.verbose==True:
			print("Periastron assumed to be at 0 degrees.")
			print(" ")

	if args.omdot and args.omdot!=0.0:
		omdot=args.omdot*(np.pi/180)*(1/365)
		fit_omdot=True
		if args.verbose==True:
			print("Periastron advance rate: {} º/yr".format(omdot))
			print("Accounting for periastron advance.")
			print(" ")
	else:
		omdot=0
		fit_omdot=False
		if args.verbose==True:
			print("Periastron advance not accounted for.")
			print(" ")

sky_position=SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
telescope_location=EarthLocation(lon=21.443888889*u.degree, lat=-30.71105556*u.degree, height=1086.6*u.m)
lon_hour=str(21.443888889*24/180)

print("")
print("event MJD, UTC range, LST range:")
write_file=open("times.txt","w")
write_file.write("MJD(event),UTC(start),UTC(end),LST(start),LST(end)\n")

k=start

while k<=end:
	if phase_def=="true" and fit_omdot==False:
		passage_mjd=get_times(phase*2*np.pi,p,e,om,t0,k)
	if phase_def=="true" and fit_omdot==True:
		passage_mjd=get_times_with_omdot(phase*2*np.pi,p,e,om,omdot,t0,k)
	if phase_def=="mean":
		passage_mjd=get_times_from_mean(phase,p,t0,k)
	passage=Time(passage_mjd,format='mjd',scale='tai')+delay
#	passage_utc=passage.utc.to_value('iso', 'date_hms')
#	passage_lst=passage.sidereal_time('mean',longitude="21.41111111")
#	local_sky_coordinates_event=sky_position.transform_to(AltAz(obstime=passage,location=telescope_location))
	local_sky_coordinates_start=sky_position.transform_to(AltAz(obstime=passage-prior*u.hour,location=telescope_location))
	local_sky_coordinates_end=sky_position.transform_to(AltAz(obstime=passage+posterior*u.hour,location=telescope_location))
	if local_sky_coordinates_start.alt>20*u.degree and local_sky_coordinates_end.alt>20*u.degree:
		passage_utc_start=(passage-prior*u.hour).utc.to_value('iso', 'date_hms')
		passage_utc_end=(passage+posterior*u.hour).utc.to_value('iso', 'date_hms')
		passage_lst_start=(passage-prior*u.hour).sidereal_time('mean',longitude=lon_hour)
		passage_lst_end=(passage+posterior*u.hour).sidereal_time('mean',longitude=lon_hour)	
		print(str(passage)+",	"+passage_utc_start+" - "+passage_utc_end+",	"+str(passage_lst_start)+" - "+str(passage_lst_end))
		write_file.write(str(passage)+","+passage_utc_start+","+passage_utc_end+","+str(passage_lst_start)+","+str(passage_lst_end)+"\n")
	k=k+1

write_file.close()