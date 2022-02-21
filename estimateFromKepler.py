import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
import argparse
import sys
#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 31/10/2021

def mass_equation(M,p_orb,a,mpulsar):
	p_orb=p_orb*24*3600
	a=299792458*a
	return M**3-(mpulsar+M)**2*4*np.pi**2*a**3/(p_orb**2*6.67408e-11)/1.9891e30

def mass_equation_derivative(M,p_orb,a,mpulsar):
	p_orb=p_orb*24*3600
	a=299792458*a
	return 3*M**2-(mpulsar+M)*8*np.pi**2*a**3/(p_orb**2*6.67408e-11)/1.9891e30

def periastron_advance(p_orb,ecc,Mtot):
	Mtot=Mtot*1.9891e30
	p_orb=p_orb*24*3600
	cnt=3*(6.67408e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(-5/3)/(1-ecc**2)
	return cnt*Mtot**(2/3)

def orbital_decay(p_orb,ecc,Mchirp):
	Mchirp=Mchirp #In kg already
	p_orb=p_orb*24*3600
	cnt=(192*np.pi/5)*(6.67408e-11/299792458**3)**(5/3)*(p_orb/(2*np.pi))**(-5/3)*(1+(73/23)*ecc**2+(37/96)*ecc**4)/(1-ecc**2)**(7/2)
	return cnt*Mchirp

def axis_decay(p_orb,ecc,Mxdot):
	p_orb=p_orb*24*3600
	Mxdot=Mxdot #In kg already
	cnt=(64/5)*299792458*(6.67408e-11/299792458**3)**2*(2*np.pi/p_orb)**2*((1+(73/23)*ecc**2+(37/96)*ecc**4)/(1-ecc**2)**(7/2))
	return cnt*Mxdot

def einstein_delay(p_orb,ecc,Mgamma):
	Mgamma=Mgamma
	p_orb=p_orb*24*3600
	gamma=ecc*(6.67408e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(1/3)
	return gamma*Mgamma

def shapiro_delay_full(Mcomp,s):
	Mcomp=Mcomp*1.9891e30
	r=6.67408e-11*Mcomp/299792458**3
	return -2*r*(np.log(1-s)-np.log(1+s))

def shapiro_delay_third(Mcomp,s):
	Mcomp=Mcomp*1.9891e30
	r=6.67408e-11*Mcomp/299792458**3
	stig=s/(1+np.sqrt(1-s**2))
	return (4*r*(stig**3)*(2/3-2*(stig**2)/5),r*(stig**3),stig)

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
		print("")
	else:
		sys.exit("The ephemeris file doesn't include all the necessary parameters.")
	return f0,p_orb,x,ecc,omega,periastron

parser=argparse.ArgumentParser(description="Take in an orbital model or orbital parameters and estimate companion mass and some post-Keplerian effects.")
parser.add_argument("--ephemeris",help="Fitorbit-format ephemeris. If given all other inputs will be ignored.")
parser.add_argument("-p","--period",type=float,help="Orbital period in days.")
parser.add_argument("-x","--axis",type=float,help="Projected semimajor axis in ls.")
parser.add_argument("-e","--eccentricity",type=float,help="Excentricity. If not given, assumed 0. Needed for higher order estimations.")
parser.add_argument("--mpulsar",type=float,help="Mass of pulsar (no uncertanties taken.). Default: 1.35 solar masses.")
parser.add_argument("--mcompanion",help="Mass of companion (mass+/-uncertainty). Use it to skip mass estimations.")
parser.add_argument("--mtotal",help="Total mass (mass+/-uncertainty). Must be used with --mcompanion. If used it has precedence over '--mpulsar'")
parser.add_argument("-i","--inclination",type=float,help="Force custom inclination angle into the mass function where applicable (degrees).")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

print(" ")

if args.ephemeris:
	(f0,p_orb,x,ecc,omega,periastron)=read_ephemeris(args.ephemeris)
else:
	if args.period:
		p_orb=args.period
		if args.verbose==True:
			print("Orbital period: {} days".format(p_orb))
			print(" ")
	else:
		sys.exit("Please specify orbital period in days with -p")
	if args.axis:
		x=args.axis
		if args.verbose==True:
			print("Projected axis: {} ls".format(x))
			print(" ")
	else:
		sys.exit("Please specify projected semimajor axis in days with -x")
	if args.eccentricity:
		ecc=args.eccentricity
		if args.verbose==True:
			print("Eccentricity: {} ls".format(x))
			print(" ")
	else:
		ecc=0
		if args.verbose==True:
			print("Eccentricity assumed to be 0")
			print(" ")

if args.mpulsar:
	mpulsar=args.mpulsar
	if args.verbose==True:
		print("Assumed pulsar mass: {} solar masses".format(mpulsar))
else:
	mpulsar=1.35
	if args.verbose==True:
		print("Assumed pulsar mass: 1.35 solar masses.")
		print(" ")

mcomp_from_massfunction=True

if args.mcompanion and args.mtotal:

	mcomp_str=args.mcompanion
	mcomp=float(mcomp_str.split("+/-")[0])
	mcomp_max=mcomp+float(mcomp_str.split("+/-")[1])
	mcomp_min=mcomp-float(mcomp_str.split("+/-")[1])

	mtot_str=args.mtotal
	mtot=float(mtot_str.split("+/-")[0])
	mtot_max=mtot+float(mtot_str.split("+/-")[1])
	mtot_min=mtot-float(mtot_str.split("+/-")[1])
	
	mpulsar=mtot-mcomp
	mpulsar_err=np.sqrt(float(mcomp_str.split("+/-")[1])**2+float(mtot_str.split("+/-")[1])**2)
	mpulsar_max=mpulsar+mpulsar_err
	mpulsar_min=mpulsar-mpulsar_err

	if args.verbose==True:
		print("Companion mass: {} solar masses".format(mcomp_str))
		print("Total mass: {} ls".format(mtotal_str))
		print("Pulsar mass overwritten to {}+/-{} solar masses".format(mpulsar,mpulsar_err))
		print(" ")

	mcomp_from_massfunction=False

elif args.mcompanion:

	mcomp_str=args.mcompanion
	mcomp=float(mcomp_str.split("+/-")[0])
	mcomp_max=mcomp+float(mcomp_str.split("+/-")[1])
	mcomp_min=mcomp+float(mcomp_str.split("+/-")[1])
	
	mpulsar_max=mpulsar
	mpulsar_min=mpulsar

	mtot=mcomp+mpulsar
	mtot_max=mcomp_max+mpulsar
	mtot_min=mcomp_min+mpulsar

	if args.verbose==True:
		print("Companion mass: {} solar masses".format(mcomp_str))
		print("Total mass: {}+/-{} ls".format(mtotal_str,mtotal_max-mtotal_min))
		print(" ")

	mcomp_from_massfunction=False

if mcomp_from_massfunction==True:

	if args.inclination:

		s=np.sin(args.inclination*np.pi/180)

		if args.verbose==True:
			print("Computing companion mass at inclination angle of {}º.".format(args.inclination))
		a_custom=x/s		
		mcomp_custom=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,a_custom,mpulsar))
		print("Companion mass at inclination {}º: {}".format(args.inclination,mcomp_custom))
		print(" ")
		mass_function_custom=4*np.pi**2*(299792458*a_custom)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30

	if args.verbose==True:
		print("Computing minimum companion mass.")
	mcomp_min=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,x,mpulsar))
	if args.inclination==True:
		print("Minimum companion mass at inclination 90º: {}".format(mcomp_min))
		print(" ")
	mass_function_min=4*np.pi**2*(299792458*x)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30

	if args.verbose==True:
		print("Computing median companion mass at inclination angle of 60º.")
	a_median=x/np.sin(60*(np.pi/180))
	mcomp_median=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,a_median,mpulsar))
	if args.inclination==True:
		print("Median companion mass at inclination 60º: {}".format(mcomp_median))
		print(" ")
	mass_function_median=4*np.pi**2*(299792458*a_median)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30

	if args.verbose==True:
		print("Computing companion mass at inclination angle of 45º.")
	a_max=x/np.sin(45*(np.pi/180))
	mcomp_max=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,a_max,mpulsar))
	if args.inclination==True:
		print("Companion mass at inclination 45º: {}".format(mcomp_max))
		print(" ")
	mass_function_max=4*np.pi**2*(299792458*a_max)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30

	if args.verbose==True:
		print("Plotting companion mass results.")
		print(" ")
	masses=np.linspace(0,3*mcomp_max/2,1000)
	masses_equation=masses**3/(mpulsar+masses)**2

	plt.plot(masses,masses_equation,"c-",label="$M_1$ = {} M$_\odot$".format(mpulsar))

	if args.inclination:
		plt.hlines(mass_function_custom,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--",label="$P_o$ = {} d, $a\\times sin(i = {}$º$)$ = {} ls".format(round(p_orb,2),args.inclination,round(x,2)))
		plt.vlines(mcomp_custom,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--",label="$M_2$ (i = {}$º$) = {} M$_\odot$".format(args.inclination,round(mcomp_custom,3)))
		plt.plot([mcomp_custom],[mass_function_custom],"ro")
	else:
		plt.hlines(mass_function_min,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--")
		plt.hlines(mass_function_median,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--",label="$P_o$ = {} d, $a\\times sin(i = 90, 60, 45$º$)$ = {} ls".format(round(p_orb,2),round(x,2)))
		plt.hlines(mass_function_max,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--")
		plt.vlines(mcomp_min,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--")
		plt.vlines(mcomp_median,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--",label="$M_2 (i = 90, 60, 45$º$)$ = {},{},{} M$_\odot$".format(round(mcomp_min,3),round(mcomp_median,3),round(mcomp_max,3)))
		plt.vlines(mcomp_max,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--")
		plt.plot([mcomp_min,mcomp_median,mcomp_max],[mass_function_min,mass_function_median,mass_function_max],"ro")

	plt.xlabel("Companion mass (M$_\odot$)")
	plt.ylabel("Mass function (M$_\odot$)")
	plt.xlim(0,3*mcomp_max/2)
	plt.ylim(np.min(masses_equation),np.max(masses_equation))
	plt.legend()
	plt.grid()
	plt.show()

	if args.inclination:
		mcomp=mcomp_custom
	else:
		mcomp=np.array([mcomp_min,mcomp_median,mcomp_max])

	mtot=mpulsar+mcomp
	if args.verbose==True:
		print("Computing periastron advance out of {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))
	omegadot=periastron_advance(p_orb,ecc,mtot)*(180/np.pi)*(24*3600*365)
	print("Expected rate of periastron advance: {} º/yr".format(omegadot))
	print(" ")

	mchirp=(mpulsar*mcomp*1.9891e30**2)/((mpulsar+mcomp)*1.9891e30)**(1/3)
	if args.verbose==True:
		print("Computing orbital decay for {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))
	pdot=orbital_decay(p_orb,ecc,mchirp)
	print("Expected rate of orbital decay: {} s/s".format(pdot))
	print(" ")

#	mxdot=(1.9891e30**2)*(mpulsar*mcomp**2)/mtot
	if args.verbose==True:
		print("Computing axis decay for {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))
#	xdot=axis_decay(p_orb,ecc,mxdot)/299792458
	xdot=x*pdot/p_orb
	print("Expected rate of axis decay: {} ls/s".format(xdot))
	print(" ")

	mgamma=(mcomp*1.9891e30)*((mpulsar+2*mcomp)*1.9891e30)/((mtot*1.9891e30)**(4/3))
	if args.verbose==True:
		print("Computing magnitude of Einstein delay for {} + {} = {} solar masses".format(mpulsar,mcomp,mtot))
	delay=einstein_delay(p_orb,ecc,mgamma)*1e6
	print("Expected magnitude of Einstein delay: {} us".format(delay))
	print(" ")

	if args.inclination:

		if args.verbose==True:
			print("Computing Saphiro delay for i = {}º".format(args.inclination))
		max_delay=shapiro_delay_full(mcomp_custom,np.sin(args.inclination*(np.pi/180)))
		(third_delay,h3,stig)=shapiro_delay_third(mcomp_custom,np.sin(args.inclination*(np.pi/180)))
		print("Expected maximum magnitude of Shapiro delay for i = {}º: {} us".format(args.inclination,max_delay*1e6))
		print("Expected magnitude of third-order Shapiro delay for i = {}º: {} us".format(args.inclination,third_delay*1e6))
		print("Expected orthometric amplitude (h3) for i = {}º: {} us".format(args.inclination,h3*1e6))
		print("Expected orthometric ratio (stig) for i = {}º: {}".format(args.inclination,stig))

	else:

		if args.verbose==True:
			print("Computing Saphiro delay for i = 85 and 75º")
		a_85=x/np.sin(85*(np.pi/180))
		mcomp_85=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,a_85,mpulsar))
		a_75=x/np.sin(75*(np.pi/180))
		mcomp_75=newton(mass_equation,mpulsar,fprime=mass_equation_derivative,args=(p_orb,a_75,mpulsar))
		max_delay_85=shapiro_delay_full(mcomp_85,np.sin(85*(np.pi/180)))
		max_delay_75=shapiro_delay_full(mcomp_85,np.sin(75*(np.pi/180)))
		(third_delay_85,h3_85,stig_85)=shapiro_delay_third(mcomp_85,np.sin(85*(np.pi/180)))
		(third_delay_75,h3_75,stig_75)=shapiro_delay_third(mcomp_75,np.sin(75*(np.pi/180)))
		print("Expected maximum magnitude of Shapiro delay for i = 85 and 75º: {} and {} us".format(max_delay_85*1e6,max_delay_75*1e6))
		print("Expected magnitude of third-order Shapiro delay for i = 85 and 75º: {} and {} us".format(third_delay_85*1e6,third_delay_75*1e6))
		print("Expected orthometric amplitude (h3) for i = 85 and 75º: {} and {} us".format(h3_85*1e6,h3_75*1e6))
		print("Expected orthometric ratio (stig) for i = 85 and 75º: {} and {}".format(stig_85,stig_75))

elif mcomp_from_massfunction==False:

	if args.verbose==True:
		print("Computing periastron advance.")
	omegadot=periastron_advance(p_orb,ecc,mtot)*(180/np.pi)*(24*3600*365)
	omegadot_max=periastron_advance(p_orb,ecc,mtot_max)*(180/np.pi)*(24*3600*365)
	omegadot_min=periastron_advance(p_orb,ecc,mtot_min)*(180/np.pi)*(24*3600*365)
	print("Estimated rate of periastron advance: {}+{}-{} º/yr".format(omegadot,omegadot_max-omegadot,omegadot-omegadot_min))
	print(" ")

	mchirp=(mpulsar*mcomp*1.9891e30**2)/((mpulsar+mcomp)*1.9891e30)**(1/3)
	mchirp_max=(mpulsar_max*mcomp_max*1.9891e30**2)/((mpulsar_max+mcomp_max)*1.9891e30)**(1/3)
	mchirp_min=(mpulsar_min*mcomp_min*1.9891e30**2)/((mpulsar_min+mcomp_min)*1.9891e30)**(1/3)
	if args.verbose==True:
		print("Computing orbital decay (period).")
	pdot=orbital_decay(p_orb,ecc,mchirp)
	pdot_max=orbital_decay(p_orb,ecc,mchirp_max)
	pdot_min=orbital_decay(p_orb,ecc,mchirp_min)
	print("Estimated rate of orbital decay: {}+{}-{} s/s".format(pdot,pdot_max-pdot,pdot-pdot_min))
	print(" ")

#	mxdot=(1.9891e30**2)*(mpulsar*mcomp**2)/mtot
#	mxdot_max=(1.9891e30**2)*(mpulsar_max*mcomp_max**2)/mtot_max
# 	mxdot_min=(1.9891e30**2)*(mpulsar_min*mcomp_min**2)/mtot_min
	if args.verbose==True:
 		print("Computing orbital decay (axis).")
#	xdot=axis_decay(p_orb,ecc,mxdot)/299792458
#	xdot_max=axis_decay(p_orb,ecc,mxdot_max)/299792458
#	xdot_min=axis_decay(p_orb,ecc,mxdot_min)/299792458
	xdot=x*pdot/p_orb
	xdot_max=x*pdot_max/p_orb
	xdot_min= x*pdot_min/p_orb
	print("Estimated rate of axis decay: {}+{}-{} ls/s".format(xdot,xdot_max-xdot,xdot-xdot_min))
	print(" ")

	mgamma=(mcomp*1.9891e30)*((mpulsar+2*mcomp)*1.9891e30)/((mtot*1.9891e30)**(4/3))
	mgamma_max=(mcomp_max*1.9891e30)*((mpulsar_max+2*mcomp_max)*1.9891e30)/((mtot_max*1.9891e30)**(4/3))
	mgamma_min=(mcomp_min*1.9891e30)*((mpulsar_min+2*mcomp_min)*1.9891e30)/((mtot_min*1.9891e30)**(4/3))	
	if args.verbose==True:
		print("Computing magnitude of Einstein delay.")
	delay=einstein_delay(p_orb,ecc,mgamma)*1e6
	delay_max=einstein_delay(p_orb,ecc,mgamma_max)*1e6
	delay_min=einstein_delay(p_orb,ecc,mgamma_min)*1e6
	print("Estimated magnitude of Einstein delay: {}+{}-{} us".format(delay,delay_max-delay,delay-delay_min))
	print(" ")

	if args.verbose==True:
		print("Computing the inclination angle.")
	s=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67408e-11)*((mtot**2)/(1.9891e30*mcomp**3)))**(1/3)
	s_min_pre=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67408e-11)*((mtot_max**2)/(1.9891e30*mcomp_max**3)))**(1/3)
	s_max_pre=299792458*x*((1/(3600*24*p_orb)**2)*(4*np.pi**2/6.67408e-11)*((mtot_min**2)/(1.9891e30*mcomp_min**3)))**(1/3)
	s_min=min(s_max_pre,s_min_pre)
	s_max=max(s_max_pre,s_min_pre)
	s_max=min(s_max,1)
	print("Estimated inclination angle: {}+{}-{} º".format(np.arcsin(s)*(180/np.pi),(np.arcsin(s_max)-np.arcsin(s))*(180/np.pi),(np.arcsin(s)-np.arcsin(s_min))*(180/np.pi)))
	print(" ")

	max_delay=shapiro_delay_full(mcomp,s)
	max_delay_max=shapiro_delay_full(mcomp_max,s_max)
	max_delay_min=shapiro_delay_full(mcomp_min,s_min)
	(third_delay,h3,stig)=shapiro_delay_third(mcomp,s)
	(third_delay_max,h3_max,stig_max)=shapiro_delay_third(mcomp_max,s_max)
	(third_delay_min,h3_min,stig_min)=shapiro_delay_third(mcomp_min,s_min)
	print("Estimated maximum magnitude of Shapiro delay: {}+{}-{} us".format(max_delay*1e6,(max_delay_max-max_delay)*1e6,(max_delay-max_delay_min)*1e6))
	print("Estimated magnitude of third-order Shapiro delay: {}+{}-{} us".format(third_delay*1e6,(third_delay_max-third_delay)*1e6,(third_delay-third_delay_min)*1e6))
	print("Estimated orthometric amplitude (h3): {}+{}-{} us".format(h3*1e6,(h3_max-h3)*1e6,(h3-h3_min)*1e6))
	print("Estimated orthometric ratio (stig): {}+{}-{}".format(stig,stig_max-stig,stig-stig_min))
