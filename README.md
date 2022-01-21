# CandyCracker

Hello! This software has been made to *crack* the secrets of binary systems discovered by the MGPS-TRAPUM pulsar survey, but it is apliable to any other pulsar survey.

## What is this repo?

This is a set of tools that allow a user to make preliminary estimates of orbital parameters out of barycentric data (epochs,periods,accelerations), plot orbital solutions and more! Useful for when you have a set of difficult points from a tricky pulsar thaat you want to solve. To add examples for the usage.

## What is NOT this repo?

This is not fitorbit! For a full implementation of fitorbit, go visit for example https://github.com/gdesvignes/pyfitorbit. However, if for example your pulsar has an orbital period smaller than the observing cadence, fitorbit will not help you too much as it needs a good period guess to start with. With ```estimateOrbit.py``` you can get this first estimate and send it to fitorbit for better luck.

## Requirements

python3, numpy, matplotlib, scipy, gatspy and astropy so far.
A tempo2 installation if you use dracula2.py.

## Contents

### estimateOrbit.py:

This snippet contains three methods to estimate the orbital period of a binary given barycentric measurements of period/accelerations. In particular:

- The ellipse method from Freire, Kramer & Lyne (2000).

- The roughness estimate from Bhattacharyya & Nityananda (2008). This particular implementation has had very significant contributions form Ewan Barr, Shalini Sengupta and David Champion. It implements an additional normalization by orbital phase during the computation of the roughness value and a check of permutations during folding to avoid skiping steps. In the future, it may also include the period errors in the analysis as well.

- A Lomb-Scargle periotrogram periodicity search.

The required packages are: numpy, matplotlib and gatspy (for the Lomb-Scargle periotogram).

### plotOrbit.py:

This one also takes barycentric data, plus a simple timing model (F0/P0 + Keplerian parameters) and plots the data along an analytical solution of the orbit in the period-phase, period-time or period-acceleration spaces. In the future, it may include a predictor mode to simulate data points in those spaces as well. Required packages: numpy and matplotlib.

### estimateFromKepler.py:

This snipped takes in Keplerian parameters from an orbital model and computes the mass function and estimates the extent of possible post-Keplerian effects. Useful to asses the nature of the system and the measurability of parameters at an early stage. Required packages: numpy, scipy and matplotlib

### dracula2.py:

I was tempted to call this snipped ```phaseConnect.py```, as it is a code to achieve a phase-connected timing solution from an initial ```tempo2``` parameter file and a set of ToAs. However, this can also be described as a ```tempo2``` implementation of the original ```Dracula``` algorithm implemented in bash for ```tempo``` (https://github.com/pfreire163/Dracula, Freire & Ridolfi 2018). Therefore, it is now called ```dracula2.py```. As such, it works very much in the very same way. Requirements: numpy and an installation of ```tempo2```. Only the basic ```tempo2``` software is called through the subprocess module, so no extra packages or plug-ins are needed.

### constrainMass.py:

This -still incomplete- snippet contrains the masses of your system given Keplerian and post-Keplerian parameters. It doesn't draw full mass and inclination diagrams, it just gives some values assuming general relativity. It only requires numpy.

### computeTimes.py:

This snippet allows you to compute observing times for MeerKAT of events happening at certain sky coordinates (e.g., periastron or conjunction at altitude higher than 20 degrees) given orbital parameters or a tempo2 ephemeris. It has some flexibility and it can account for omdot. Eventually, more telescopes sites can be added. It requires numpy, scipy and (of course) astropy.

## To do's of this repo:

To implement errors into the roughness and ellipse estimators in ```estimateOrbit.py```. To implement a Bayesian method for the ellipse. To put the functions into a callable library instead of the headers of each script.

## Questions, inquiries and feedback

Just write to mcbernadich@mpifr-bonn.mpg.de if you have any.