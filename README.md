# Local Group Halo Properties

Inferring the total mass, and other properties, of the halos of the most massive local group galaxies, MW, M31 and M33.

We take the observational constraints on the sky position of, distance to, radial velocity and proper motion of M31, M33, and the LMC, and transform them into constraints on the M31-centric distances, radial velocities and tangential velocities. We do this my drawing samples from the posterior PDFs for the observables, transforming them to the M31 frame, and then approximating the resulting multi-variate PDF with a Gaussian Mixture Model (GMM).

We then compute the value of this PDF for each quadruplet of Consuelo N-body simulation halos, which were pre-selected to be local group-like in their isolation and very approximate mass: this likelihood is used to weight each quad, as we build up a posterior PDF for the halo properties (principally masses and concentrations) based on the Consuelo-sampled, GMM-approximated prior. We also compare our local group mass with some simple timing argument estimates. 

This work builds on [our (Busha et al's) earlier work on the Milky Way halo](http://adsabs.harvard.edu/abs/2011ApJ...743...40B), and on [Gonzalez,  Kravtsov and Gnedin's initial investigation of the Local Group system](http://adsabs.harvard.edu/abs/2014ApJ...793...91G).

### Products

* **Paper: [Williamson et al in prep.](https://www.overleaf.com/read/sxjdrsbqfwqt)**
* [Poster](https://www.overleaf.com/read/bkwpxstxpqns)
* [Williamson 2015, Stanford Senior Undergraduate Thesis](https://www.overleaf.com/read/pxsmhcmkwdmh)

### Contact

* Marc Williamson (KIPAC, NYU)
* Phil Marshall (KIPAC)
* Yao-Yuan Mao (KIPAC, U Pitt)
* Risa Wechsler (KIPAC)
* Michael Busha

This is research in progress, and the code is available for re-use under GPL v2. If you make use of any of the ideas or code in this repository in your own research we ask you to cite us, and acknowledge this code repository, please. The current reference is "Williamson 2015, Williamson et al, in prep." Willamson 2015 is Marc Williamson's undergrad senior thesis, which can be obtained at [this link](https://www.overleaf.com/read/pxsmhcmkwdmh). If you have comments or questions about this work, or would like to get involved, [please do get in touch via the issues!](https://github.com/drphilmarshall/LocalGroupHaloProps/issues)


### Notebooks, demos etc

* [Research notes](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/notes.ipynb)
* Notebook: [Triplets: MW, M31, M33](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/triplets.ipynb)
* Notebook: [Triplets, with the GMM prior](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/gmm_prior.ipynb)
* Notebook: [Pairs: MW, M31, with the GMM prior](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/gmm_pair_prior.ipynb)
* Notebook: [Pairs, straight from Consuelo](http://nbviewer.ipython.org/github/drphilmarshall/LocalGroupHaloProps/blob/master/pairs.ipynb)


### Using the code

Add the module to your `$PYTHONPATH` and your `$PATH`, and the scripts should run. You'll need the Consuelo halo data of course - we have uploaded a small sample for test purposes.

You will also need the following library modules:

    import numpy,sys,getopt,string,subprocess

and also `astropy`, and Yao-Yuan Mao's "helpers" module:

    git clone git@bitbucket.org:yymao/helpers.git

