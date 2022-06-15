# rns
Rapidly Rotating Neutron Star (RNS)

Rapidly Rotating Neutron Star
=============================

RNS is a code written by [Nikolaos Stergioulas](
http://www.astro.auth.gr/~niksterg/) which constructs models of rapidly rotating, relativistic, compact stars using tabulated equations of state which are supplied by the user. Please direct questions about this program to either [Sharon Morsink](mailto:morsink@phys.ualberta.ca) or [Nikolaos Stergioulas.](mailto:niksterg@aei-potsdam.mpg.de)

The code was mainly developed while Stergioulas and Morsink worked at UWM, Nikolaos as a graduate student and Sharon as a postdoc. Currently Nikolaos Stergioulas is an assistant professor at the Aristotle University of Thessaloniki (Greece) and Sharon Morsink is an Associate Professor at University of Alberta (Canada).

Description
-----------

The code is based on the method developed by [Komatsu, Eriguchi & Hachisu (1989)](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1989MNRAS.237..355K) and modifications introduced by [Cook, Shapiro & Teukolsky (1994).](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1994ApJ...422..227C) It can compute individual models as well as sequences of fixed mass, rest mass, angular velocity or angular momentum models. All models assume uniform rotation. You can read more about this code in [Stergioulas and Friedman (1995)](http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=1995ApJ...444..306S) and in [Nozawa, Stergioulas, Gourgoulhon & Eriguchi (1998)](http://xxx.lanl.gov/abs/gr-qc/9804048) and references therein.

A [list of publications](rns/publications.html) which use the rns code is available.

Sources
-------

*   Version 1 This version of the code is best for users who simply wish to input an equation of state and a parameter (such as angular velocity, mass, etc.) and have the computer output a list of the star's properties.
    
    *   [rns.v1.1c.tar.gz](rns/source/rns.v1.1c.tar.gz) Nick's 1996 version.
    *   [rns.v1.1d.tar.gz](rns/source/rns.v1.1d.tar.gz) Includes modifications by Sharon Morsink (1997): Prints the star's quadrupole moment (Based on method by [Laarakkers and Poisson (1997)](http://xxx.lanl.gov/abs/gr-qc/9709033)); Units for polytropic stars are corrected.
    
      
    
*   Version 2 This version of the code is more suitable if you wish to use the results (such as the metric) in another application. The code is modularized so that it is fairly easy to call up from other programs. An example program is included. (Please email [Sharon Morsink](mailto:morsink@phys.ualberta.ca) with suggestions which will make this more user friendly.)
    
    *   [rns.v2.0.tar.gz](rns/source/rns.v2.0.tar.gz) Apr. 13, 1999.
    
      
    
*   A [tarred gzip file](rns/source/eos.tar.gz) of nuclear equation of state files are available. An index to the [available EOS](rns/source/eos/EOS.INDEX) files is included.  
      
    
*   The code, makefile, equation of state files and user's manual are available for [download.](rns/source/)

HTML User Manual
----------------

*   [Version 1.1x](rns/manual1/)
*   [Version 2.0](rns/manual2/)

Links
-----

*   [Nikolaos Stergioulas' web page, with more information and videos about Rotating Relativistic Stars](http://www.astro.auth.gr/~niksterg/research.html)
*   [Sharon Morsink's web page](http://fermi.phys.ualberta.ca/~morsink/)

Thanks .....
------------

*   Greg Cook - for supplying some equations of state files
*   Shai Ayal - for pointing out an error in a previous version of 2.0
