import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt
import argparse
import sys
#Writen by Miquel Colom Bernadich i la mare que el va parir. Last update: 31/10/2021

def mass_equation(M,p_orb,x,mpulsar):
	p_orb=p_orb*24*3600
	x=299792458*x
	return M**3-(mpulsar+M)**2*4*np.pi**2*x**3/(p_orb**2*6.67408e-11)/1.9891e30

def mass_equation_derivative(M,p_orb,x,mpulsar):
	p_orb=p_orb*24*3600
	x=299792458*x
	return 3*M**2-(mpulsar+M)*8*np.pi**2*x**3/(p_orb**2*6.67408e-11)/1.9891e30

def periastron_advance(p_orb,ecc,Mtot):
	Mtot=Mtot*1.9891e30
	p_orb=p_orb*24*3600
	cnt=3*(6.67408e-11/299792458**3)**(2/3)*(p_orb/(2*np.pi))**(-5/3)/(1-ecc**2)
	return cnt*Mtot**(2/3)

def orbital_decay(p_orb,ecc,Mchirp):
	Mchirp=Mchirp
	p_orb=p_orb*24*3600
	cnt=-(192*np.pi/5)*(6.67408e-11/299792458**3)**(5/3)*(p_orb/(2*np.pi))**(-5/3)*(1+(73/23)*ecc**2+(37/96)*ecc**4)/(1-ecc**2)**(7/2)
	return cnt*Mchirp

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

parser=argparse.ArgumentParser(description="Take in an orbiral model or orbital parameters and constrain companion mass and some post-Keplerian effects.")
parser.add_argument("--ephemeris",help="Fitorbit-format ephemeris. If given all other inputs will be ignored.")
parser.add_argument("-p","--period",type=float,help="Orbital period in days")
parser.add_argument("-a","--axis",type=float,help="Projected semimajor axis in ls.")
parser.add_argument("-e","--eccentricity",type=float,help="Excentricity. If not given, assumed 0. Needed for higher order estimations.")
parser.add_argument("--mpulsar",type=float,help="Mass of pulsar. Default: 1.4 solar masses")
parser.add_argument("--mcompanion",type=float,help="Mass of companion. Use it to skip mass estimations")
parser.add_argument("-v","--verbose",action="store_true")
args = parser.parse_args()

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
		sys.exit("Please specify projected semimajor axis in days with -a")
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
		print("Pulsar mass: {} ls".format(mpulsar))
else:
	mpulsar=1.4
	if args.verbose==True:
		print("Pulsar mass assumed to be 1.4 solar masses.")
		print(" ")

if args.mcompanion:
	mcomp=args.mcompanion
	if args.verbose==True:
		print("Companion mass: {} solar masses".format(mcomp))
		print(" ")
else:

	if args.verbose==True:
		print("Computing minimum companion mass.")
	mcomp_min=newton(mass_equation,1.4,fprime=mass_equation_derivative,args=(p_orb,x,mpulsar))
	print("Minimum companion mass at inclination 90º: {}".format(mcomp_min))
	print(" ")

	if args.verbose==True:
		print("Computing median companion mass at inclination angle of 60º.")
	a_median=x/np.sin(60*(np.pi/180))
	mcomp_median=newton(mass_equation,1.4,fprime=mass_equation_derivative,args=(p_orb,a_median,mpulsar))
	print("Median companion mass at inclination 60º: {}".format(mcomp_median))
	print(" ")

	if args.verbose==True:
		print("Computing companion mass at inclination angle of 45º.")
	a_max=x/np.sin(45*(np.pi/180))
	mcomp_max=newton(mass_equation,1.4,fprime=mass_equation_derivative,args=(p_orb,a_max,mpulsar))
	print("Companion mass at inclination 45º: {}".format(mcomp_max))
	print(" ")

	if args.verbose==True:
		print("Plotting companion mass results.")
		print(" ")
	masses=np.linspace(0,3*mcomp_max/2,1000)
	masses_equation=masses**3/(mpulsar+masses)**2
	mass_function_min=4*np.pi**2*(299792458*x)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30
	mass_function_median=4*np.pi**2*(299792458*a_median)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30
	mass_function_max=4*np.pi**2*(299792458*a_max)**3/((24*3600*p_orb)**2*6.67408e-11)/1.9891e30

	plt.plot(masses,masses_equation,"c-",label="$M_1$ = {} M$_\odot$".format(mpulsar))
	plt.hlines(mass_function_min,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--")
	plt.hlines(mass_function_median,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="-",label="$P_o$ = {} d, $a\\times sin(i = 90, 60, 45º)$ = {} ls".format(round(p_orb,2),round(x,2),60))
	plt.hlines(mass_function_max,max(mcomp_min-2,0),mcomp_max+2,color="blue",linestyles="--")
	plt.vlines(mcomp_min,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--")
	plt.vlines(mcomp_median,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="-",label="$M_2$ (i = 90, 60, 45º) = {},{},{} M$_\odot$".format(round(mcomp_min,2),round(mcomp_median,2),round(mcomp_max,2)))
	plt.vlines(mcomp_max,np.min(masses_equation),np.max(masses_equation),color="red",linestyles="--")
	plt.plot([mcomp_min,mcomp_median,mcomp_max],[mass_function_min,mass_function_median,mass_function_max],"ro")

	plt.xlabel("Trial mass (M$_\odot$)")
	plt.ylabel("Mass function (M$_\odot$)")
	plt.xlim(0,3*mcomp_max/2)
	plt.ylim(np.min(masses_equation),np.max(masses_equation))
	plt.legend()
	plt.grid()
	plt.show()

	mcomp=np.array([mcomp_min,mcomp_median,mcomp_max])

mtot=mpulsar+mcomp
if args.verbose==True:
	print("Computing periastron advance out of {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))
omegadot=periastron_advance(p_orb,ecc,mtot)*(180/np.pi)*(24*3600*365)
print("Estimated rate of periastron advance: {} º/yr".format(omegadot))
print(" ")
if args.verbose==True:
	print("Computing periastron advance out of {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))

mchirp=(mpulsar*mcomp*1.9891e30**2)/((mpulsar+mcomp)*1.9891e30)**(1/3)
if args.verbose==True:
	print("Computing orbital decay for {} + {} = {} solar masses.".format(mpulsar,mcomp,mtot))
pdot=orbital_decay(p_orb,ecc,mchirp)*(24*3600*365)
pdot_tempo2=orbital_decay(p_orb,ecc,mchirp)/(24*3600)*1e12
print("Estimated rate of periastron advance: {}e-12 d/s = {} s/yr".format(pdot_tempo2,pdot))
print(" ")


