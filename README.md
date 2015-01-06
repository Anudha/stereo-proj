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

## Requirements

* Python 2.7 with Matplotlib, Numpy, Tkinter, PIL (available through [Enthought Canopy](https://store.enthought.com/downloads/) or [Anaconda](http://continuum.io/downloads) distributions for instance).
* Mac users need to install [Active Tcl](http://www.activestate.com/activetcl/downloads).
* Run the script python /my-path-to-the-script/stereo-proj.py

# User guide

## Interface

[img1](/img1.png)

##Plotting procedure
### From a set of diffraction vectors

* Enter the crystal structure or import it from the menu bar. The structure can be modified/added by modifying the structure.txt file. The format is: name a b c alpha beta gamma space-group. 
* Enter a diffraction vector, the tilt angle and the inclination angle determined by [Index](https://github.com/mompiou/index) for instance and press plot button below.
* Example for Al: diffraction 2: (111), tilt: 30, inclination: 147.5 (tilt axis is vertical here). Rotate here of 47° around the diffraction vector (choose a rotation increment) to bring the diffraction vector 1 in the right position. The Euler angles are given. Running the mouse over the projection gives you the inclination and tilt angles.
[img2](/img2.png)

### From the Euler angles
* Enter the angles in the "Euler angles" field and press the plot button right below.

##Additional features
* Draw the crystal directions by ticking the uvw button.
* Add a pole/direction by entering the indices in the "Add a pole" field and by clicking the "Add" button.
* Draw a plane by entering the indices in the "Add a pole" field and by clicking the "Plane" button.
* Draw equivalent poles/directions or planes by entering the indices in the "Add a pole" field and by clicking the "Symmetry" or "Sym planes" buttons, respectively.
* The "minus" buttons allow you to erase planes/pole/directions entered in the "Add a pole" field.
* Calculate interplanar spacing by entering the indices in the "Add a pole" field and by clicking the "dhkl" button.
* Draw the iso-contours of Schmid factor for a given plane by entering the indices in the "Add a pole" field and by clicking the "iso-schmid" button (straining axis along the tilt axis).
[img3](/img3.png)
* Right click on the projection to draw the closest pole/direction (maximum indice 8)
* Change the number of pole/direction drawn by increasing or decreasing the value of "d" (choose an increment)
* Draw the apparent width variation with the tilt angle for an given plane by filling the "width" field and by clicking on the "width" button.
[img4](/img4.png)
* Compute the angle between two poles by filling the angle field and by clicking on the "angle" button
* Calculate the Schmid factor by filling the "b" (Burgers vector) and "n" (slip plane) field and by clicking on the "schmid factor" button
* Make rotation along x, y, z by entering a rotation increment ant by clicking on "+/-" buttons in the "x,y,z rotation field
* For hexagonal structure, tick the "hexa" button to draw the projection with 4 indices.
* Change the layout by selecting "square/circle", color by ticking "green, "blue" or "red" and marker size by changing the value (40 default) before plotting.

## Export
* Save the projection from the menu bar (jpeg default, works with pdf)
