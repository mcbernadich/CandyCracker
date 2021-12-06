import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 31/10/2021

def fit_ellipse(periods,derivatives):

	print("Assuming circular orbit and guessing...")
	print("")

	coeff=np.polynomial.polynomial.polyfit(periods,derivatives**2,2)

	period=-coeff[1]/(2*coeff[2])
	p_orb=2*np.pi*299792458/(period*np.sqrt(-coeff[2]))
	axis=p_orb*np.sqrt(period**2-coeff[0]/coeff[2])/(2*np.pi*period)
	acc_amp=np.square(2*np.pi/p_orb)*axis*299792458
	p_amp=2*np.pi*period*axis/p_orb

	print("Spin period= {} ms".format(period*1000))
	print("Orbital period= {} days".format(p_orb/24/3600))
	print("Projected axis= {} ls".format(axis))
	phase=np.arctan(-history[2]/(acc_amp*(history[1]-period)*history[1]))
	per_epoch=history[0]-phase*p_orb/(2*np.pi*24*3600)
	print("First epoch of periastron= {} MJD".format(per_epoch[0]))
	print("")

	return (period,p_orb,axis,acc_amp,p_amp,coeff)

def fit_folding(epochs,periods,errors,period_range="none"):

	import time

	sorted_epochs=np.sort(epochs)
	interval=sorted_epochs[-1]-sorted_epochs[0]
	if period_range != "none":
		trial_porb=float(period_range.split(":")[0])
		last_porb=float(period_range.split(":")[1])
	else:
		trial_porb=0.041
		last_porb=interval
	smallest_interval=np.min(sorted_epochs[1:]-sorted_epochs[:-1])
	epochs=epochs-sorted_epochs[0]
	sorting_old=np.argsort(epochs)

	print("Searching for best orbital period in between {} and {} days...".format(trial_porb,last_porb))
	print("")

	roughness=[]
	trial_porbs=[]
	unsolvable_knots=0
	correction_attempts=0
	i=0

	start=time.time()

	# Start Kernel.

	while trial_porb <= last_porb:

		if np.mod(i,10000)==0:
			print("Current period: {} days. Current step: {} days".format(trial_porb,1e-2*(trial_porb**2)/(2*np.pi*interval)))
		folded_epochs=epochs-(epochs//trial_porb)*trial_porb
		sorting=np.argsort(folded_epochs)

		# Consistency check.
		if (sorting!=sorting_old).any() and i!=0:
			iterations=0
			good_to_go=False
			while good_to_go==False and iterations<=30:
				likelihood=(sorting==sorting_old)
				# First, check if a single alone can explain the shift (first point always the same, therefore excluded).
				if (np.roll(sorting[1:],1)==sorting_old[1:]).all():
					good_to_go=True					
				# If not, check if a single permutation can explain the difference.
				elif likelihood[np.where(likelihood==False)].size==2:
					good_to_go=True
				# If none of the previous, start going back and forth with the steps until one of the preious holds.
				if (sorting!=sorting_old).any() and good_to_go==False:
					step=step/2
					trial_porb=trial_porb-step
					folded_epochs=epochs-(epochs//trial_porb)*trial_porb
					sorting=np.argsort(folded_epochs)
				if (sorting==sorting_old).all() and good_to_go==False:
					step=step/2
					trial_porb=trial_porb+step
					folded_epochs=epochs-(epochs//trial_porb)*trial_porb
					sorting=np.argsort(folded_epochs)
				iterations=iterations+1
			if iterations>1 and good_to_go==True:
				correction_attempts=correction_attempts+1
			if iterations==31 and good_to_go==False:
				step=old_step
				trial_porb=trial_porb_old+step
				folded_epochs=epochs-(epochs//trial_porb)*trial_porb
				sorting=np.argsort(folded_epochs)
				correction_attempts=correction_attempts+1
				unsolvable_knots=unsolvable_knots+1
		# -----------------------------------------------------------------------
		
		folded_periods=periods[sorting]
		folded_phases=folded_epochs[sorting]/trial_porb
		phase_difference=np.roll(folded_phases,-1)-folded_phases
		phase_difference[ phase_difference <= 0 ] = phase_difference[ phase_difference <= 0 ] + 1
		roughness.append(np.sum((((np.roll(folded_periods,-1)-folded_periods)/phase_difference)**2)))
		trial_porbs.append(trial_porb)
		# For consistency check
		trial_porb_old=trial_porb
		sorting_old=sorting
		step=1e-2*(trial_porb**2)/(2*np.pi*interval)
		old_step=step
		# -------------------------------------------------------------------------
		trial_porb=trial_porb+step
		i=i+1

	# End Kernel.

	print("Run time: {} s".format(time.time()-start))

	roughness=np.array(roughness)
	trial_porbs=np.array(trial_porbs)
	print("Total number of trials: {}".format(np.size(roughness)))
	print("Number of correction attempts: {}".format(correction_attempts))
	print("Number of unsolvable knots: {}".format(unsolvable_knots))
	print("Succesful corrections: {}".format(correction_attempts-unsolvable_knots))

	best_porb=trial_porbs[np.argmin(roughness)]

	print("Best orbital period= {} days".format(best_porb))
	print("")

	return (best_porb,trial_porbs,roughness)

def fit_LombScargle(epochs,periods,errors,period_range="none"):

	from gatspy.periodic import LombScargle

	sorted_epochs=np.sort(epochs)
	interval=sorted_epochs[-1]-sorted_epochs[0]
	if period_range != "none":
		min_porb=float(period_range.split(":")[0])
		last_porb=float(period_range.split(":")[1])
	else:
		min_porb=0.041
		last_porb=interval
	smallest_interval=np.min(sorted_epochs[1:]-sorted_epochs[:-1])

	ls=LombScargle().fit(epochs,periods,errors)
	(trial_porbs,power)=ls.periodogram_auto(nyquist_factor=500)
	print(trial_porbs,power)
	ls.optimizer.period_range = (min_porb, last_porb)
	best_porb=ls.best_period

	return (best_porb,trial_porbs,power)

def read_ephemeris(file):
	ephemeris=open(file,"r")
	ecc=0
	omega=0
	for line in ephemeris:
		line=line.replace(" ","")
		line=line.replace("	","")
		if line[:2]=="F0":
			f0=float(line[2:])
		elif line[:2]=="PB":
			p_orb=float(line[2:])
		elif line[:2]=="A1":
			x=float(line[2:])
		elif line[:3]=="ECC":
			ecc=float(line[3:])
		elif line[:2]=="OM":
			omega=float(line[2:])
		elif line[:2]=="T0":
			periastron=float(line[2:])
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

parser=argparse.ArgumentParser(description="Take in barycentric data at each epoch and fit it to some orbital parameters.")
parser.add_argument("data",help="File with columns of 1) MJD, 2) period and 3) uncertainty (or 3) derivatives, 4) derivative uncertainty).")
parser.add_argument("-p","--period",help="Units of period. Default: 'ms'.",choices=["s","ms","s-1"])
parser.add_argument("-d","--derivative",help="Units of derivative. Default: 's/s'.",choices=["s/s","s-2","m/s2"])
parser.add_argument("-r","--range",help="Range of data lines to be read 'min:max'. Default: '0:inf'.")
parser.add_argument("-m","--method",help="Method to fit the data with",choices=["ellipse","roughness","Lomb-Scargle"])
parser.add_argument("--period_range",help="Custom orbital period search range in days for -m 'roughness' and 'Lomb-Scargle'. 'min:max' in days")
parser.add_argument("-M","--methods",help="Print details about methods.",action="store_true")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

if args.methods:
	print("Currently, there are 2 models t fit the data with:")
	print("")
	print(" - Ellipse:")
	print("   It reads the 2nd (Period), 3rd (Derivative) and 4th (Derivative uncertainty) columns of data and apply Freire, Kramer & Lyne (2000). It fits a parabola on the absolute values of acceleration in function of baricentric period. The coefficients tell you about the orbit. Very good if you have good period and derivative values and the orbit is not excentric.")
	print("")
	print(" - Roughness:")
	print("   It reads the 1st (MJD) and 2nd (Period) columns of data and apply Bhattacharyya & Nityananda (2008). It folds your data points at consecutive trial periods in between 1 hour and the duration of the data set, estimating a roughness value that should be minimal at the true period. It can take a long time to run, but is is very useful if the orbit is very eccentric or you don't have good derivative measureements.")
	print("   If '--period_range' is specified, then the orbital period search range is restricted in days to the custom values.")
	print(" - Lomb-Scargle:")
	print("   It reads the 1st (MJD) and 2nd (Period) columns of data and computes a Lomb-Scargle diagrame to search the best period in between 1 hour and the duration of the data set. Useful if you don't have good derivative measurements, much faster than the roughness algorithm, but less sensible to highly eccentric orbits.")
	print("   If '--period_range' is specified, then the orbital period search range is restricted in days to the custom values.")
	print("")
	exit()

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

if args.range:
	start=int(args.args.range.split(":")[0])
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
	if args.period_range:
		period_range=args.period_range
	if args.verbose:
		print("The data will be fitted with the following method: "+args.method)
		if model=="roughness" and args.period_range:
			print("The folding range will be of "+period_range+" days.")
		print("Write '-M' to know more about the models.")
		print("")
else:
	fit=False
	model="none"
	if args.verbose:
		print("The data will just be show in a period - epoch diagram")
		print("Give '-m' to estimate orbit with a model.")
		print("Write '-M' to know more about the models.")
		print("")

if p_units=="s-1":
	history[1]=1/history[1]
	if model!="ellipse":
		history[2]=history[2]/np.square(history[1])

if p_units=="ms":
	history[1]=history[1]/1000
	if model!="ellipse":
		history[2]=history[2]/1000

if pdot_units=="s/s" and model=="ellipse":
	history[2]=299792458*history[2]/history[1]
	history[3]=299792458*history[3]/history[1]

if pdot_units=="s-2" and model=="ellipse":
	if args.verbose:
		print("Loaded frequency derivative values:")
		print(history[2])
	history[2]=-299792458*history[2]*history[1]
	history[3]=-299792458*history[3]*history[1]
	if args.verbose:
		print("Translated acceleration values:")
		print(history[2])
		print("")

if fit==True:
	if model == "ellipse":
		(period,p_orb,axis,acc_amp,p_amp,coeff)=fit_ellipse(history[1],history[2])
		plot_fit=True
	elif model == "roughness":
		if args.period_range:
			(p_orb,trial_porbs,roughness)=fit_folding(history[0],history[1],history[2],period_range=period_range)
		else:
			(p_orb,trial_porbs,roughness)=fit_folding(history[0],history[1],history[2])
		plot_fit=True
	elif model == "Lomb-Scargle":
		if args.period_range:
			(p_orb,trial_porbs,ls_power)=fit_LombScargle(history[0],history[1],history[2],period_range=period_range)
		else:
			(p_orb,trial_porbs,ls_power)=fit_LombScargle(history[0],history[1],history[2])
	else:
		sys.exit("Please state a valid orbit estimation algorithm.")
else:
	plt.errorbar(history[0]-int(history[0,0]),history[1]*1000,history[2]*1000,fmt="o")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("epoch - {} (MJD)".format(int(history[0,0])))
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.show()

if fit == True and model == "ellipse":

	plt.errorbar(1000*history[1],history[2]**2,yerr=2*history[2]*history[3],fmt="o")
	periods=np.linspace(period-p_amp,period+p_amp,1000)
	plt.plot(1000*periods,coeff[0]+coeff[1]*periods+coeff[2]*periods**2,"c-")
	plt.plot([],[]," ",label="$P_{pulsar}=$ "+str(round(period*1000,4))+" ms,\n$P_{orb}=$ "+str(round(p_orb/(24*3600),3))+" d,\n$x=$ "+str(round(axis,3))+" ls")
	plt.ylabel("$acc^2$ (m/s$^2$)$^2$")
	plt.xlabel("$P_{bary}$ (ms)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()

	plt.errorbar(1000*history[1],history[2],yerr=history[3],fmt="o")
	phases=np.linspace(0,2*np.pi,1000)
	periods=period+p_amp*np.cos(phases)
	accs=-acc_amp*np.sin(phases)
	plt.plot(1000*periods,accs,"c-")
	plt.plot([],[]," ",label="$P_{pulsar}=$ "+str(round(period*1000,4))+" ms,\n$P_{orb}=$ "+str(round(p_orb/(24*3600),3))+" d,\n$x=$ "+str(round(axis,3))+" ls")
	plt.ylabel("$acc$ (m/s$^2$)")
	plt.xlabel("$P_{bary}$ (ms)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()

elif fit == True and model == "roughness":
	
	plt.plot(trial_porbs,np.max(roughness)/roughness)
	plt.plot([],[]," ",label="$P_{orb}=$ "+str(round(p_orb,5))+" d")
	plt.ylabel("Reciprocal $Roughness$")
	plt.xlabel("Trial orbital period (days)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.legend()
	plt.tight_layout()
	plt.show()

	epochs=history[0,:]-np.min(history[0,:])
	folded_epochs=epochs-(epochs//p_orb)*p_orb
	plt.errorbar(folded_epochs,history[1]*1000,history[2]*1000,fmt="o")
	plt.plot([],[]," ",label="$P_{orb}=$ "+str(round(p_orb,5))+" d")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("orbital phase (MJD)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()

elif fit == True and model == "Lomb-Scargle":

	plt.plot(trial_porbs,ls_power)
	plt.plot([],[]," ",label="$P_{orb}=$ "+str(round(p_orb,5))+" d")
	plt.ylabel("Lomb-Scargle power")
	plt.xlabel("Trial orbital period (days)")
	if args.period_range:
		plt.xlim(float(period_range.split(":")[0]),float(period_range.split(":")[1]))
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.legend()
	plt.tight_layout()
	plt.show()

	epochs=history[0,:]-np.min(history[0,:])
	folded_epochs=epochs-(epochs//p_orb)*p_orb
	plt.errorbar(folded_epochs,history[1]*1000,history[2]*1000,fmt="o")
	plt.plot([],[]," ",label="$P_{orb}=$ "+str(round(p_orb,5))+" d")
	plt.ylabel("$P_{bary}$ (ms)")
	plt.xlabel("orbital phase (MJD)")
	plt.title(args.data.split(".")[0].split("/")[-1])
	plt.tight_layout()
	plt.legend()
	plt.show()