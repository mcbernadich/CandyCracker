import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 31/10/2021

def frequency_fromE(E,Edot,f_0,fxom,fxomecc):
	return f_0+fxom*np.sin(E)*Edot-fxomecc*np.cos(E)*Edot #This can be vectorial too.

def frequency_derivative_fromE(E,Edot,Edotdot,fxom,fxomecc):
	return fxom*(np.sin(E)*Edotdot+np.cos(E)*np.square(Edot))-fxomecc*(np.cos(E)*Edotdot-np.sin(E)*np.square(Edot)) #This one can be vectorial too.

def period_accel(f_0,Porb,x,ecc,omega,n_points):
	#Compute eccentric anomaly and time derivatives.
	E=np.linspace(0,2*np.pi,n_points)
	pi2orb=2*np.pi/(Porb*86400)
	Edot=pi2orb/(1-ecc*np.cos(E))
	Edotdot=-(ecc*np.sin(E)*Edot**2)/(1-ecc*np.cos(E))
	#Compute constants for the frequency equations.
	fxom=f_0*x*np.sin(np.deg2rad(omega))
	fxomecc=f_0*x*np.cos(np.deg2rad(omega))*np.sqrt(1-ecc**2)
	#Compute periods and accelerations.
	periods=1/frequency_fromE(E,Edot,f_0,fxom,fxomecc)
	accelerations=-299792458*periods*frequency_derivative_fromE(E,Edot,Edotdot,fxom,fxomecc)
	return periods,accelerations #This can be vectorial too.

def time_period(f_0,Porb,x,ecc,omega,periastron,n_points):
	#Compute eccentric anomaly and time derivatives.
	E=np.linspace(0,2*np.pi,n_points)
	pi2orb=2*np.pi/(Porb*86400)
	Edot=pi2orb/(1-ecc*np.cos(E))
	Edotdot=-(ecc*np.sin(E)*Edot**2)/(1-ecc*np.cos(E))
	#Compute time.
	time=(E-ecc*np.sin(E))/pi2orb+periastron*86400
	#Compute constants for the frequency equations.
	fxom=f_0*x*np.sin(np.deg2rad(omega))
	fxomecc=f_0*x*np.cos(np.deg2rad(omega))*np.sqrt(1-ecc**2)
	#Compute periods.
	periods=1/frequency_fromE(E,Edot,f_0,fxom,fxomecc)
	return time/86400,periods #This can be vectorial too.

def time_der(f_0,Porb,x,ecc,omega,periastron,n_points):
	#Compute eccentric anomaly and time derivatives.
	E=np.linspace(0,2*np.pi,n_points)
	pi2orb=2*np.pi/(Porb*86400)
	Edot=pi2orb/(1-ecc*np.cos(E))
	Edotdot=-(ecc*np.sin(E)*Edot**2)/(1-ecc*np.cos(E))
	#Compute time.
	time=(E-ecc*np.sin(E))/pi2orb+periastron*86400
	#Compute constants for the frequency equations.
	fxom=f_0*x*np.sin(np.deg2rad(omega))
	fxomecc=f_0*x*np.cos(np.deg2rad(omega))*np.sqrt(1-ecc**2)
	#Compute periods.
	periods=1/frequency_fromE(E,Edot,f_0,fxom,fxomecc)
	derivatives=-frequency_derivative_fromE(E,Edot,Edotdot,fxom,fxomecc)*periods**2
	return time/86400,derivatives #This can be vectorial too.

def time_acc(f_0,Porb,x,ecc,omega,periastron,n_points):
	#Compute eccentric anomaly and time derivatives.
	E=np.linspace(0,2*np.pi,n_points)
	pi2orb=2*np.pi/(Porb*86400)
	Edot=pi2orb/(1-ecc*np.cos(E))
	Edotdot=-(ecc*np.sin(E)*Edot**2)/(1-ecc*np.cos(E))
	#Compute time.
	time=(E-ecc*np.sin(E))/pi2orb+periastron*86400
	#Compute constants for the frequency equations.
	fxom=f_0*x*np.sin(np.deg2rad(omega))
	fxomecc=f_0*x*np.cos(np.deg2rad(omega))*np.sqrt(1-ecc**2)
	#Compute periods.
	periods=1/frequency_fromE(E,Edot,f_0,fxom,fxomecc)
	accelerations=-299792458*periods*frequency_derivative_fromE(E,Edot,Edotdot,fxom,fxomecc)
	return time/86400,accelerations #This can be vectorial too.

def read_ephemeris(file):
	ephemeris=open(file,"r")
	ecc=0
	omega=0
	for line in ephemeris:
		line = line.strip().split()
		if line[0]=="F0":
			f0=float(line[1])
		elif line[0]=="PB":
			p_orb=float(line[1])
		elif line[0]=="A1":
			x=float(line[1])
		elif line[0]=="ECC":
			ecc=float(line[1])
		elif line[0]=="OM":
			omega=float(line[1])
		elif line[0]=="T0":
			periastron=float(line[1])
	if f0 and p_orb and x and periastron:
		print("Pulsar parameters loaded from {}:".format(file))
		print("- Pulsar frequency: {} Hz".format(f0))
		print("- Orbital period: {} days".format(p_orb))
		print("- Projected axis: {} ls".format(x))
		print("- Eccentricity: {}".format(ecc))
		print("- Periastron angle: {} degrees".format(omega))
		print("- Periastron passage: {} (MJD)".format(periastron))
	else:
		sys.exit("The ephemeris file doesn't include all the necessary parameters.")
	return f0,p_orb,x,ecc,omega,periastron

parser=argparse.ArgumentParser(description="Take in barycentric data at each epoch, an orbital model, and plot for your eyes to enjoy.")
parser.add_argument("--data",help="File with columns of 1) MJD, 2) period and 3) uncertainty (or 3) derivatives, 4) derivative uncertainty).")
parser.add_argument("--ephemeris",help="Fitorbit-format ephemeris.")
parser.add_argument("-p","--period",help="Units of period values. Default: 'ms'.",choices=["s","ms","s-1"])
parser.add_argument("-d","--derivative",help="Units of derivative values. Default: 's/s'.",choices=["s/s","s-2","m/s2"])
parser.add_argument("-r","--range",help="Range of data lines to be read 'min:max'. Default: '0:inf'.")
parser.add_argument("-m","--method",help="Method to plot the data with.",choices=["p-acc","p-phase","acc-phase","d-phase","pulse_drift"])
parser.add_argument("--mjd",type=float,help="Plot a value point at that MJD")
#parser.add_argument("--residuals",action="store_true",help="Print and store the period residuals of data vs model")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

if args.period:
	p_units=args.period
else:
	p_units="ms"
if args.verbose==True:
	print("Period units are "+p_units+".")
	print("")

if args.derivative:
	pdot_units=args.derivative
else:
	pdot_units="s/s"
if args.verbose==True:
	print("Derivative units are "+pdot_units+".")
	print("")

if args.data:

	if args.range:
		start=int(args.range.split(":")[0])
		end=int(args.range.split(":")[1])
		history=np.loadtxt(data).T[:,start:end+1]
		if args.verbose==True:
			print("Reading ephemeris from {} to {} from file {}".format(start,end,args.data))
			print("")
	else:
		history=np.loadtxt(args.data).T
		if args.verbose==True:
			print("Reading all ephemeris from file {}".format(args.data))
			print("")

if args.method:
	fit=True
	model=args.method
	if args.verbose:
		print("The data will be plotted with a model in the following display: "+args.method)
		print("")
else:
	fit=False
	model="none"
	if args.verbose:
		print("The data will just be show in a period - epoch diagram")
		print("Give '-m' to plot orbit with an orbital model.")
		print("")

if args.data:

	if p_units=="s-1":
		history[1]=1/history[1]
		if model!="p-acc":
			history[2]=history[2]/np.square(history[1])

	if p_units=="ms":
		history[1]=history[1]/1000
		if model!="p-acc":
			history[2]=history[2]/1000

	if pdot_units=="s/s" and model=="p-acc":
		history[2]=299792458*history[2]/history[1]
		history[3]=299792458*history[3]/history[1]

	if pdot_units=="s-2" and model=="p-acc":
		if args.verbose:
			print("Loaded frequency derivative values:")
			print(history[2])
		history[2]=-299792458*history[2]*history[1]
		history[3]=-299792458*history[3]*history[1]
		if args.verbose:
			print("Translated acceleration values:")
			print(history[2])
			print("")

(f0,p_orb,x,ecc,omega,periastron)=read_ephemeris(args.ephemeris)

if fit == False and args.data:

	plt.errorbar(history[0]-int(history[0,0]),history[1]*1000,history[2]*1000,fmt="o")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("epoch - {} (MJD)".format(int(history[0,0])))
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.show()

elif fit == True and model == "p-acc":

	(periods,accelerations)=period_accel(f0,p_orb,x,ecc,omega,2000)

	if args.data:
		plt.errorbar(1000*history[1],history[2],yerr=history[3],fmt="o")
	plt.plot(1000*periods,accelerations,"c-")
	plt.plot([],[]," ",label="$P_0$ = {} ms,\n$P_B$ = {} d, $x$ = {} ls,\n$ecc$ = {}, $\omega$ = {}º".format(round(1000/f0,3),round(p_orb,3),round(x,3),round(ecc,3),round(omega,3)))
	plt.ylabel("$acc$ (m/s$^2$)")
	plt.xlabel("$P_{bary}$ (ms)")
#	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.legend()
	plt.tight_layout()
	plt.show()

elif fit == True and model == "p-phase":

	if args.data:
		first_epoch=np.min(history[0,:])
		epochs=history[0,:]-first_epoch
		reduced_periastron=periastron-first_epoch
		folded_epochs=epochs-(epochs//p_orb)*p_orb
		reduced_periastron=reduced_periastron-(reduced_periastron//p_orb)*p_orb
		folded_epochs[ folded_epochs < reduced_periastron ] = folded_epochs[ folded_epochs < reduced_periastron ] + p_orb
		plt.errorbar((folded_epochs-reduced_periastron)/p_orb,history[1]*1000,history[2]*1000,fmt="o")
	else:
		reduced_periastron=periastron

	(time,periods)=time_period(f0,p_orb,x,ecc,omega,reduced_periastron,2000)
	plt.plot((time-reduced_periastron)/p_orb,1000*periods,"c-")
	plt.plot([],[]," ",label="$P_0$ = {} ms, $P_B$ = {} d,\n$x$ = {} ls, $ecc$ = {},\n$\omega$ = {}º, $T_0$= {} MJD".format(round(1000/f0,3),round(p_orb,3),round(x,3),round(ecc,3),round(omega,3),round(periastron,3)))
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("Mean anomaly (phase)")
#	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.legend()
	plt.tight_layout()
	plt.xlim(0,1)
	plt.show()

elif fit == True and (model == "d-phase" or model == "acc-phase"):

	if args.data:
		first_epoch=np.min(history[0,:])
		epochs=history[0,:]-first_epoch
		reduced_periastron=periastron-first_epoch
		folded_epochs=epochs-(epochs//p_orb)*p_orb
		reduced_periastron=reduced_periastron-(reduced_periastron//p_orb)*p_orb
		folded_epochs[ folded_epochs < reduced_periastron ] = folded_epochs[ folded_epochs < reduced_periastron ] + p_orb
		plt.errorbar((folded_epochs-reduced_periastron)/p_orb,history[1]*1000,history[2]*1000,fmt="o")
	else:
		reduced_periastron=periastron

	if model == "d-phase":
		(time,derivatives)=time_der(f0,p_orb,x,ecc,omega,reduced_periastron,2000)
		plt.ylabel("P$_1$ (s/s)")
	if model == "acc-phase":
		(time,derivatives)=time_acc(f0,p_orb,x,ecc,omega,reduced_periastron,2000)
		plt.ylabel("$a$ (m/$s^2$)")

	plt.plot((time-reduced_periastron)/p_orb,derivatives,"c-")
	plt.plot([],[]," ",label="$P_0$ = {} ms, $P_B$ = {} d,\n$x$ = {} ls, $ecc$ = {},\n$\omega$ = {}º, $T_0$= {} MJD".format(round(1000/f0,3),round(p_orb,3),round(x,3),round(ecc,3),round(omega,3),round(periastron,3)))

	if args.mjd:
		reduced_mjd=p_orb*((args.mjd-periastron)/p_orb-(args.mjd-periastron)//p_orb)+periastron
		reduced_time=abs(time-reduced_mjd)
		closest_time=time[np.argpartition(reduced_time,0)[0]]
		closest_2nd_time=time[np.argpartition(reduced_time,1)[1]]
		closest_derivative=derivatives[np.argpartition(reduced_time,0)[0]]
		closest_2nd_derivative=derivatives[np.argpartition(reduced_time,1)[1]]
		if min(closest_time,closest_2nd_time)==closest_2nd_time:
			derivative_mjd=closest_2nd_derivative+(reduced_mjd-closest_2nd_time)*(closest_derivative-closest_2nd_derivative)/(closest_time-closest_2nd_time)
		if min(closest_time,closest_2nd_time)==closest_time:
			derivative_mjd=closest_derivative+(reduced_mjd-closest_time)*(closest_2nd_derivative-closest_derivative)/(closest_2nd_time-closest_time)
		print("")
		print("MJD= "+str(args.mjd)+" days lays at orbital phase "+str((reduced_mjd-periastron)/p_orb)+", with a(P1)="+str(derivative_mjd)+" m/s2(s/s)")
		plt.plot((reduced_mjd-periastron)/p_orb,derivative_mjd,"ro")

	plt.xlabel("Mean anomaly (phase)")
#	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.legend()
	plt.tight_layout()
	plt.xlim(0,1)
	plt.show()