## University of Luxembourg
# Non-Parametric Methods in Econometrics (2019)
Supplementary files for prof. Gautam Tripathi’s course ‘Non-parametric methods‘ at the University of Luxembourg maintained by Andreï V. Kostyrka.

R files in this project:
* _smoothing-functions.R_ — basic functions to work with kernels;
* _smoothing-simulation-01-univariate.R_ — how to use these functions for one-dimensional smoothing;
* _smoothing-simulation-02-multivariate.R_ — how to use these functions for multi-dimensional smoothing;
* _smoothing-simulation-03-distribution.R_ — examining theoretical properties of kernel estimators and their distributions.

PDF files in this project (with TeX sources):
* _typeset-code.pdf_ — code with syntax highlighting and images generate by it;
* _handouts.pdf_ — formulæ for kernel convolutions.

You will need ImageMagick installed on your computer if you want to generate animated plots from simulation 03. If you do not succeed, the GIF animations are in the `/output/` folder.

This code was tested on Linux (Mint 19.1 and Arch) and Windows 10. Since grid search and simulations can be run in parallel, Linux or Mac are highly advocated since vanilla R on Windows does not have parallel support.
