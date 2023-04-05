# Lévy Nonnegative Matrix Factorization

This repository contains the code related to the Lévy NMF algorithm for robust nonnegative source separation. It was published in two papers (see [below](#references)) that you're encouraged to check and to cite if you use the related content.


## Using this code

### Setup

Even though this code was initially developed with Matlab, we've adapted it to Octave. To fully use it, you need several Octave packages, which we can install as follows:

	sudo apt install liboctave-dev
	sudo apt install octave-control
	sudo apt install octave-signal
	sudo apt install octave-statistics

Note that you can ignore `statistics` if you don't intend to plot the results. This repository also provides and uses the robust PCA code. Please cite the [corresponding paper](https://posenhuang.github.io/papers/RPCA_Separation_ICASSP2012.pdf) if you use it.

## Data

We provide in the `data/` folder:

- The data used for the experiment on fluorescence spectroscopy, which has been provided by Cyril Gobinet from [Université de Reims](https://www.univ-reims.fr/biospect/presentation/membres/membres,22804,37993.html).
- Short guitar excerpts from the [IDMT-SMT-GUITAR](https://zenodo.org/record/7544110) database. Please cite the corresponding paper if you use it.

In addition, our experiments use the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de/) for music separation. Download it, and unzip it in the `data` folder (or change the dataset path accordingly).


### Reproducing the papers' results

Run the files in the `scripts/` folder to reproduce the various experiments from our papers (fitting impulsive noise, fluorescence spectroscopy, music inpainting, and music accompaniment enhancement).

Several additional scripts are provided in the `analysis/` folder to explore the behavior of Lévy NMF (robustness to initialization, convergence of the algorithms, majorize-equalization performance). 


## References

The results from this project have been published in the following papers:

- P. Magron, R. Badeau, A. Liutkus, [Lévy NMF for robust nonnegative source separation](https://hal.science/hal-01548488), Proc. IEEE WASPAA, 2017.
- P. Magron, R. Badeau, A. Liutkus, [Lévy NMF : un modèle robuste de séparation de sources non-négatives](https://hal.science/hal-01540484), Actes du XXVIe Colloque GRETSI, 2017 (in French).


