LocalGroupHaloProps
===================

Inferring the total mass, and other properties, of the halos of the most massive local group galaxies, MW, M31 and M33.

The Demo shows the most polished, up to date results of the project. It is meant to organize the plots that will eventually be used in the paper. The Notes go into more detail on bug fixes and code exploration. It includes experimenting with sk-learn methods, figuring out the triangle.py package, etc.

Demo [here](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/demo.ipynb)

Notes [here](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/notes.ipynb)

### Program

We take the observational constraints on the sky position of, distance to, radial velocity and proper motion of M31 and M33, and transform them into constraints on the M31-centric distances, radial velocities and tangential velocities of M33 and the Milky Way. We do this my drawing samples from the posterior PDFs for the observables, transforming them to the M31 frame, and then approximating the resulting 6-D PDF.

We then compute the value of this PDF at each triplet of Consuelo halos, which were pre-selected to be local group-like in their isolation and very approximate mass: this likelihood is ised to weight each triplet, as we build up a posterior PDF for the halo properties (principally masses and shapes) based on the Consuelo-sampled prior.

We also compare our local group mass with some simple timing argument estimates. This work builds on [our (Busha et al's) earlier work on the Milky Way halo](http://adsabs.harvard.edu/abs/2011ApJ...743...40B), and on [Gonzalez,  Kravtsov and Gnedin's initial investigation of the Local Group system](http://adsabs.harvard.edu/abs/2014ApJ...793...91G).

### Getting started

Add the module to your `$PYTHONPATH` and your `$PATH`, and the scripts should run. You'll need the Consuelo halo data of course - we have uploaded a small sample for test purposes.

You will also need the following library modules:

    import numpy,sys,getopt,string,subprocess

and also `astropy`, and Yao-Yuan Mao's "helpers" module:

    git clone git@bitbucket.org:yymao/helpers.git


### Authors

* Marc Williamson (KIPAC)
* Phil Marshall (KIPAC)
* Yao-Yuan Mao (KIPAC)
* Risa Wechsler (KIPAC)
* Michael Busha (KIPAC)


### License, Citation etc

This code is available for re-use under GPL v2. However, if you make use of it in your research we ask you to cite our forthcoming paper and acknowledge this code repository, please. The paper reference is "Williamson et al, in preparation" - the work documented here will form Marc Williamson's undergrad thesis project. You can follow the development of this thesis at [overleaf](https://www.overleaf.com/read/pxsmhcmkwdmh).
