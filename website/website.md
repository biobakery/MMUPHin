---
* Install MMUPHin in R with `devtools::install_bitbucket("biobakery/mmuphin@master")`
* Or download the [package source code](https://bitbucket.org/biobakery/mmuphin/get/HEAD.zip)
* For detailed installation and usage help, please consult the [MMUPHin manual](https://bitbucket.org/biobakery/waafle/src/master/README.md)
---

MMUPHin is an R package for meta-analysis tasks of microbiome cohorts. It has function interfaces for a) covariate-controlled batch- and cohort effect adjustment, b) meta-analysis differential abundance testing, c) meta-analysis unsupervised discrete structure (clustering) discovery, and d) meta-analysis unsupervised continuous structure discovery. 

## Citing MMUPHin

A manuscript describing MMUPHin, with its application in a collection of inflammatory bowel disease studies is currently in prep:

|Siyuan Ma, Dmitry Shungin, Himel Mallick, Melanie Schirmer, Long Nguyen, Raivo Kolde, Eric Franzosa, Hera Vlamakis, Ramnik Xavier, Curtis Huttenhower *The landscape of novel lateral gene transfer events in the human microbiome.*|
|---|

In the meantime, if you use WAAFLE or the datasets provided below in your work, please cite this website: http://huttenhower.sph.harvard.edu/mmuphin.


## Quick-start guide

### Requirements

MMUPHin is an R package and its interfaces are provided as R functions. It requires the following dependencies. Note that all packages except for MaAsLin2 will be installed automatically the first time you invoke the installation command. For installation of MaAsLin2, please refer to [its webpage](http://huttenhower.sph.harvard.edu/maaslin2).

* [Base R 3.5.0](https://www.r-project.org/) or newer
* CRAN packages: `metafor`,`fpc`,`igraph`,`ggplot2`,`cowplot`,`ggrepel`
* [MaAsLin2](http://huttenhower.sph.harvard.edu/maaslin2)

### Installation
* If not already available, install `devtools`:
	* `> install.packages("devtools")`
* Once MaAsLin2 was installed, run the following command to install MMUPHin (all the other dependencies will be installed automatically for you):
	* `> devtools::install_bitbucket("biobakery/mmuphin@master")`

### Performing meta-analysis with MMUPHin

* Detailed tutorial to come. For the time being, refer to the help manuals of MMUPHin's functions for their usage (e.g., run `?MMUPHin::adjust.batch` in R)