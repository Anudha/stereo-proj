Stereo-proj
===========

Stereo-proj is a python script software that plots the stereographic projection of a crystal given its lattice parameters (a,b,c,α,β,γ) for either a set of 3 Euler angles (ϕ1,Φ,ϕ2) or a pole and an angle of rotation. You can additionnally:
* Add a specific pole (h,k,l) /direction [u,v,w] or equivalent ones (special 4 indices mode is available for hexagonal structures)
* Draw a specific plane or a set of equivalent planes
* Click and find poles/directions
* Select the number of poles/directions to draw
* Rotate the crystal around x,y and z directions
* Calculate interplanar spacing dhkl, angles between directions/planes
* Plot iso-Schmid factor lines for a given plane, and calculate Schmid factor
* Draw the apparent variation of a plane width with tilt angle
* Save plot and import structure menu

Requirements
=============
* Require python 2.7 with Matplotlib, Numpy, Tkinter, PIL (use Enthought Canopy or Anaconda).
* For Mac users, you need to install Active Tcl.
* Run the script python /path/stereo.py

