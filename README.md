LocalGroupHaloProps
===================

Inferring the mass etc of the local group galaxy halos (MW, M31 and M33).

Demo [here](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/demo.ipynb)

### Program

We take the observational constraints on the sky position of, distance to, radial velocity and proper motion of M31 and M33, and transform them into constraints on the M31-centric distances, radial velocities and tangential velocities of M33 and the Milky Way. We do this my drawing samples from the posterior PDFs for the observables, transforming them to the M31 frame, and then approximating the resulting 6-D PDF.

We then compute the value of this PDF at each triplet of Consuelo halos, which were pre-selected to be local group-like in their isolation and very approximate mass: this likelihood is ised to weight each triplet, as we build up a posterior PDF for the halo properties (principally masses and shapes) based on the Consuelo-sampled prior.

We also compare our local group mass with some simple timing argument estimates.

### Getting started

Add the module to your `$PYTHONPATH` and your `$PATH`, and the scripts should run. You'll need the Consuelo halo data of course - we'll try and upload a small sample for test purposes?

You will also need the following library modules:

    import numpy,atpy,sys,getopt,string,subprocess,pyfits

and also `astropy`...

### Authors

* Phil Marshall (KIPAC)
* Michael Busha (KIPAC)
* Marc Williamson (KIPAC)


### License, Citation etc

This code is available for re-use under GPL v2. However, if you make use of it in your research we ask you to cite our forthcoming paper and acknowledge this code repository, please. The paper reference is "Williamson et al, in preparation" - the work documented here will form Marc Williamson's undergrad thesis project.
