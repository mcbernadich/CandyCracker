import numpy as np
import subprocess as sub
from os.path import exists
import glob
import argparse

def read_chi2r(tempo2logs):
	for line in tempo2logs.split("\n"):
		if line[0:11]=="Fit Chisq =":
			chi2r=float(line.split("=")[3].split("	")[0])
		if line[0:25]=="Number of fit parameters:":
			nfit=int(line.split(":")[1])
		if line[0:25]=="Number of points in fit =":
			npoints=int(line.split("=")[1])
	return chi2r,npoints-nfit

def change_inclination(parFile,M2,SINI):

	r=6.67408e-11*1.9891e30*M2/299792458**3
	stig=SINI/(1+np.sqrt(1-SINI**2))
	h3=r*(stig**3)

	par_read=open(parFile,"r")
	parFile_new=parFile.split(".")[0]+"_"+str(SINI)+".par"
	par_write=open(parFile_new,"w")
	for line in par_read:
		chunks = line.strip().split()
		if chunks[0]=="M2":
			par_write.write("M2 "+str(M2)+"\n")
		elif chunks[0]=="H3":
			par_write.write("H3 "+str(h3)+"\n")
		elif chunks[0]=="STIG":
			par_write.write("STIG "+str(stig)+"\n")
		else:
			par_write.write(line)

	return parFile_new

parser=argparse.ArgumentParser(description="Take in a DDGR or DDH file and make a chi2r map on fixed inclination values.")
parser.add_argument("-p","--parameter",help="Tempo2 parameter file. It doesn't matter if there is M2 or h3 and stig values, they will all be changed.")
parser.add_argument("-t","--tim",help="Tempo2 tim file.")
parser.add_argument("-i","--inclinations",help="List of values of SINI, M2 (or CHI2).")
parser.add_argument("-m","--model",help="Select the tempo2 model.",choices=["DDH","DDGR"])
parser.add_argument("--chi2_list",help="File with the list of SINI values and CHI2 values")
args = parser.parse_args()

if args.parameter:
	parFile=args.parameter
if args.tim:
	timFile=args.tim
if args.inclinations:
	incls=np.loadtxt(args.inclinations,delimiter=",")
sub.run(["mkdir","chi2_parFiles"])

i=0
if args.model=="DDH" or args.model=="DDGR":

	chi2_list=[]
	chi2_file=open("chi2_list"+args.model+".txt","w")
	
	for row in incls:

		SINI=row[0]
		M2=row[1]
		parFile_new=change_inclination(parFile,M2,SINI)

		sub.run(["tempo2","-f",parFile_new,timFile,"-outpar",parFile_new],stdout=sub.DEVNULL)
		sub.run(["tempo2","-f",parFile_new,timFile,"-outpar",parFile_new],stdout=sub.DEVNULL)
		sub.run(["tempo2","-f",parFile_new,timFile,"-outpar",parFile_new],stdout=sub.DEVNULL)
		a=sub.run(["tempo2","-f",parFile_new,timFile,"-outpar",parFile_new],stdout=sub.PIPE)

		(chi2r,nfree)=read_chi2r(a.stdout.decode())
		print(SINI,chi2r*nfree)
		chi2_list.append(chi2r*nfree)
		sub.run(["mv",parFile_new,"chi2_parFiles"],stdout=sub.DEVNULL)

		i=i+1

		if args.model=="DDH":
			r=6.67408e-11*1.9891e30*row[1]/299792458**3
			stig=row[0]/(1+np.sqrt(1-row[0]**2))
			h3=r*(stig**3)
			chi2_file.write("{},{},{},{}\n".format(row[0],h3,stig,chi2r*nfree))
		if args.model=="DDGR":
			chi2_file.write("{},{}\n".format(row[0],row[1],chi2r*nfree))

	chi2_file.close()
	chi2_list=np.array(chi2_list)

else:

	incls=np.loadtxt(args.chi2_list,delimiter=",").T
	chi2_list=incls[1]
	sini=incls[0]

best_arg=np.argmin(chi2_list)

best_sini=sini[best_arg]

reduced_chi2_list=chi2_list-np.min(chi2_list)-1
upper_sini=sini[np.argmin(np.abs(reduced_chi2_list[0:best_arg]))]
lower_sini=sini[np.argmin(np.abs(reduced_chi2_list[best_arg:]))]

print("Inclination angle range from chi2 mapping:")
print("")
print("sin(i)= "+str(best_sini)+" +"+str(upper_sini-best_sini)+" -"+str(best_sini-lower_sini))
print("i= "+str(np.arcsin(best_sini)*(180/np.pi))+" +"+str((np.arcsin(upper_sini)-np.arcsin(best_sini))*(180/np.pi))+" -"+str((np.arcsin(best_sini)-np.arcsin(lower_sini))*(180/np.pi))+" degrees")
print("")