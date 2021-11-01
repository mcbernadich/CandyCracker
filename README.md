# CandyCracker

Hello! This software has been made to *crack* the secrets of binary systems discovered by the MGPS-TRAPUM pulsar survey, but it is apliable to any other pulsar survey.

## What is this repo?

This is a set of tools that allow a user to make preliminary estimates of orbital parameters out of barycentric data (epochs,periods,accelerations), plot orbital solutions and more! Useful for when you have a set of difficult points from a tricky pulsar thaat you want to solve.

## What is NOT this repo?

This is not fitorbit! For a full implementation of fitorbit, go visit for example https://github.com/gdesvignes/pyfitorbit. However, if for example your pulsar has an orbital period smaller than the observing cadence, fitorbit will not help you too much as it need a good period guess to start with. With ```estimateOrbit.py``` you can get this first estimate and send it to fitorbit for better luck.

## Requirements

python3, numpy, matplotlib, scipy and gatspy so far.

## To do's of this repo:

To implement errors into the roughness and ellipse estimators in ```estimateOrbit.py```. To implement a Bayesian method for the ellipse. To put the functions into a callable library instead of the headers of each script.

## Questions, inquiries and feedback

Just write to mcbernadich@mpifr-bonn.mpg.de if you have any.