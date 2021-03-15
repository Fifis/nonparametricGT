## University of Luxembourg
# nonparametricGT: Simple Functions for the Non-Parametric Methods in Econometrics course (2021)
Supplementary files for prof. Gautam Tripathi’s Ph.D.-level course ‘Non-parametric methods‘ at the University of Luxembourg written and maintained by Andreï V. Kostyrka.

This is most likely to be the final version because I am using more general versions of these functions in my Ph.D. thesis, so these simple functions should be used for didactic purposes. This small package will help you answer these questions:

- How can I get a simple, no-nonsense, vanilla Parzen—Rosenblatt density estimator?
- How can I perform non-parametric smoothing, no-nonsense, vanilla Parzen—Rosenblatt density estimator?

## Installation

No compilation is required! On Windows, Mac, and Linux, just run these lines:

    install.packages("devtools") # If you do not have it already
    library(devtools)
    install_github("Fifis/nonparametricGT", subdir = "package")
    library(nonparametricGT)

It will pull the _nonparametricGT_ package from this repository (subdirectory _package_) and install it on your computer as a local package.

In case you would like to examine all the functions manually, see

R files in this project:
* _package/R/smoothing-functions.R_ — basic functions to work with kernels;
* _smoothing-simulation-01-univariate.R_ — how to use these functions for one-dimensional smoothing;
* _smoothing-simulation-02-multivariate.R_ — how to use these functions for multi-dimensional smoothing;
* _smoothing-simulation-03-distribution.R_ — examining theoretical properties of kernel estimators and their distributions.
* _smoothing-simulation-04-mixed.R_ — non-parametric density estimation and regression when the data types are mixed
* _smoothing-simulation-05-efficient.R_ — implementation of P. M. Robinson’s (1987, Econometrica) semi-parametrically efficient estimator

PDF files in this project (with TeX sources):
* _handouts.pdf_ — formulæ for kernel convolutions.
* _nonparametricGT_0.1.1.pdf_ — package help in PDF format

You will need ImageMagick installed on your computer if you want to generate animated plots from simulation 03. If you do not succeed, the GIF animations are in the `/output/` folder. Besides that, to produce high-quality MP4 videos of the 3D shapes, you need to have `ffmpeg` installed, as well as a HEVC (H.265) encoder (`sudo apt install x265`, `sudo pacman -S x265` etc.).

This code was tested on Linux (Mint 20.1 and Arch) and Windows 10 (except for the video rendering part). Since grid search and simulations can be run in parallel, Linux or Mac are highly advocated since vanilla R on Windows does not have parallel support.
