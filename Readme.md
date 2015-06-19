[![DOI](https://zenodo.org/badge/9647/chumo/Gating_with_adatoms.svg)](http://dx.doi.org/10.5281/zenodo.14781)

Gating_with_adatoms
===================

This is a program to visualize the electrostatic potential 
induced by charged defects on surfaces, as they are placed by the user.

For a detailed explanation on the physics behind this, see the articles:

*Vertical manipulation of native adatoms on the InAs(111)A surface*
Yang, J.; Nacci, C.; Martínez-Blanco, J.; Kanisawa, K. & Fölsch, S.
[Journal of Physics: Condensed Matter, 2012, 24, 354008](http://dx.doi.org/10.1088/0953-8984/24/35/354008)

*Emergent multistability in assembled nanostructures*
Yang, J.; Erwin, S.; Kanisawa, K.; Nacci, C. & Fölsch, S.
[Nano Letters, 2011, 11, 2486](http://dx.doi.org/10.1021/nl2009444)

This repo contains a Python version and an HTML version:

### Python
The graphical user interface is a Matplotlib figure which is able to handle events
like mouse clicks and mouse movements.

It has been tested in Python 2.7.6 and Python 3.4.2.

For the moment, the program deals only with the case of adatoms adsorbed on 
the vacancy sites of the 2x2 reconstruction of the InAs(111)A semiconductor surface. 

### HTML
This is a web-app where the interactivity has been implemented using the javascript visualization library d3.js.

In this version it is possible to define a lattice of adsorption sites (two vectors, each with modulus and angle).

Apart from [d3.js](http://d3js.org/), a copy of two javascript libraries are included and used, namely [PNGlib.js](http://www.xarg.org/2010/03/generate-client-side-png-files-using-javascript/) and [Bessel.JS](https://github.com/SheetJS/bessel). A guided tour through the app is also implemented, for which [jQuery.js](https://jquery.com/) and [Bootstrap-tour.js](http://bootstraptour.com/) are used.

You can try the web-app via [rawgit](https://rawgit.com/chumo/Gating_with_adatoms/master/WebApp/Gating_with_adatoms.html).

Enjoy!

Jesús Martínez-Blanco
