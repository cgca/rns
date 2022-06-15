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

Links
-----

*   [Nikolaos Stergioulas' web page, with more information and videos about Rotating Relativistic Stars](http://www.astro.auth.gr/~niksterg/research.html)
*   [Sharon Morsink's web page](http://fermi.phys.ualberta.ca/~morsink/)

Thanks .....
------------

*   Greg Cook - for supplying some equations of state files
*   Shai Ayal - for pointing out an error in a previous version of 2.0

RNS Publication List
====================

### 1\. [astro-ph/9411032](http://xxx.lanl.gov/abs/astro-ph/9411032) :

Title: COMPARING MODELS OF RAPIDLY ROTATING RELATIVISTIC STARS CONSTRUCTED BY TWO NUMERICAL METHODS  
Authors: [Nikolaos Stergioulas](http://xxx.lanl.gov/find/astro-ph/1/Nikolaos+Stergioulas/0/1/0/all/1/1), [John L. Friedman](http://xxx.lanl.gov/find/astro-ph/1/John+L%2E+Friedman/0/1/0/all/1/1)  
Comments: 12 pages, AASTEX + 10 postscript figures and postscript version of tex file, accepted for publication in the Astrophysical Journal  
Journal-ref: ApJ 444 (1995) 306

### 2\. [N. Stergioulas, PhD Thesis](ftp://pauli.phys.uwm.edu/pub/rns/nst_thesis.ps.gz) :

Title: The Structure and Stability of Rotating Relativistic Stars  
Authors: [Nikolaos Stergioulas](http://xxx.lanl.gov/find/astro-ph/1/Nikolaos+Stergioulas/0/1/0/all/1/1)  
Comments: Ph.D Thesis, UWM, Department of Physics, 1996

### 3\. [astro-ph/9608179](http://xxx.lanl.gov/abs/astro-ph/9608179) :

Title: Upper Limit Set by Causality on the Rotation and Mass of Uniformly Rotating Relativistic Stars  
Authors: [Scott Koranda](http://xxx.lanl.gov/find/astro-ph/1/Scott+Koranda/0/1/0/all/1/1) (Case Western Reserve University), [Nikolaos Stergioulas](http://xxx.lanl.gov/find/astro-ph/1/Nikolaos+Stergioulas/0/1/0/all/1/1), [John L. Friedman](http://xxx.lanl.gov/find/astro-ph/1/John+L%2E+Friedman/0/1/0/all/1/1) (University of Wisconsin- Milwaukee)  
Comments: 28 pages, LaTeX2e, 8 Postscript figures, submitted to ApJ  
Journal-ref: ApJ 488 (1997) 799

### 4\. [gr-qc/9705056](http://xxx.lanl.gov/abs/gr-qc/9705056) :

Title: Nonaxisymmetric Neutral Modes in Rotating Relativistic Stars  
Authors: [Nikolaos Stergioulas](http://xxx.lanl.gov/find/gr-qc/1/Nikolaos+Stergioulas/0/1/0/all/1/1), [John L. Friedman](http://xxx.lanl.gov/find/gr-qc/1/John+L%2E+Friedman/0/1/0/all/1/1)  
Comments: 57 pages, 9 figures  
Journal-ref: ApJ 492 (1998) 301

### 5\. [gr-qc/9709033](http://xxx.lanl.gov/abs/gr-qc/9709033) :

Title: Quadrupole moments of rotating neutron stars  
Authors: [William G. Laarakkers](http://xxx.lanl.gov/find/gr-qc/1/William+G%2E+Laarakkers/0/1/0/all/1/1), [Eric Poisson](http://xxx.lanl.gov/find/gr-qc/1/Eric+Poisson/0/1/0/all/1/1)  
Comments: ReVTeX, 7 pages, 5 figures, additional material, and references added  
Journal-ref:

### 6\. [astro-ph/9712120](http://xxx.lanl.gov/abs/astro-ph/9712120) :

Title: A General Relativistic Calculation of Boundary Layer and Disk Luminosity for Accreting Non-Magnetic Neutron Stars in Rapid Rotation  
Authors: [Arun V. Thampan](http://xxx.lanl.gov/find/astro-ph/1/Arun+V%2E+Thampan/0/1/0/all/1/1) (1), [Bhaskar Datta](http://xxx.lanl.gov/find/astro-ph/1/Bhaskar+Datta/0/1/0/all/1/1) (1 and 2) ((1)Indian Institute of Astrophysics, Bangalore, India, (2)Raman Research Institute, Bangalore, India)  
Comments: LaTeX file using mn.sty, 11 pages (including 7 figs, 1 table), uses psbox.tex. Accepted for publication in MNRAS  
Journal-ref:

### 7\. [gr-qc/9804048](http://xxx.lanl.gov/abs/gr-qc/9804048) :

Title: Construction of Highly Accurate Models of Rotating Neutron Stars - Comparison of Three Different Numerical Schemes  
Authors: [Tetsuo Nozawa](http://xxx.lanl.gov/find/gr-qc/1/Tetsuo+Nozawa/0/1/0/all/1/1), [Nikolaos Stergioulas](http://xxx.lanl.gov/find/gr-qc/1/Nikolaos+Stergioulas/0/1/0/all/1/1), [Eric Gourgoulhon](http://xxx.lanl.gov/find/gr-qc/1/Eric+Gourgoulhon/0/1/0/all/1/1), [Yoshiharu Eriguchi](http://xxx.lanl.gov/find/gr-qc/1/Yoshiharu+Eriguchi/0/1/0/all/1/1)  
Comments: 24 pages with 16 figures, submitted to Astronomy and Astrophysics  
Journal-ref: A&AS 132 (1998) 431

### 8\. [gr-qc/9805012](http://xxx.lanl.gov/abs/gr-qc/9805012) :

Title: Rotating Stars in Relativity  
Authors: [Nikolaos Stergioulas](http://xxx.lanl.gov/find/gr-qc/1/Nikolaos+Stergioulas/0/1/0/all/1/1)  
Comments: 23 pages, review article for Living Reviews in Relativity, [this http URL](http://www.livingreviews.org) replaced with final version (corrected minor typing errors, updated several sections and added new references)  
Journal-ref: Living Reviews in Relativity, Vol. 1, No. 1998-8, (1998)  

### 9\. [gr-qc/9806008](http://xxx.lanl.gov/abs/gr-qc/9806008) :

Title: Quasi-normal modes of rotating relativistic stars - neutral modes for realistic equations of state  
Authors: [Sharon M. Morsink](http://xxx.lanl.gov/find/gr-qc/1/Sharon+M%2E+Morsink/0/1/0/all/1/1), [Nikolaos Stergioulas](http://xxx.lanl.gov/find/gr-qc/1/Nikolaos+Stergioulas/0/1/0/all/1/1), [Steve R. Blattnig](http://xxx.lanl.gov/find/gr-qc/1/Steve+R%2E+Blattnig/0/1/0/all/1/1)  
Comments: 12 pages, 5 figures, submitted to ApJ  
Journal-ref:

### 10\. [astro-ph/9808227](http://xxx.lanl.gov/abs/astro-ph/9808227) :

Title: Relativistic precession around rotating neutron stars: Effects due to frame-dragging and stellar oblateness  
Authors: [Sharon M. Morsink](http://xxx.lanl.gov/find/astro-ph/1/Sharon+M%2E+Morsink/0/1/0/all/1/1), [Luigi Stella](http://xxx.lanl.gov/find/astro-ph/1/Luigi+Stella/0/1/0/all/1/1)  
Comments: 35 pages including 10 figures and 6 tables. To appear in the Astrophysical Journal  
Journal-ref:

To add or update references, email: [Sharon Morsink](mailto:morsink@phys.ualberta.ca)
