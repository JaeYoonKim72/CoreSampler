# CoreSampler: core samples/subset selection from large genotype datasets

CoreSampler is a modified version of GenoCore, which improves on the limitations of GenoCore where the user is unable to specify the number of samples to extract.

The CoreSampler extracts core samples as many numbers as user wanted, and provides a VCF file of core samples and a coverage file recorded how well the core samples are representative to the entire sample.

This program is based on the GenoCore python version availiable at https://github.com/JaeYoonKim72/GenoCore_Python.

Source code was written by Jae-Yoon Kim using Python language and supported on windows and linux platform.


# CoreSampler algorithm

The algorithm schematic of "CoreSet", the main submodule, is as follows.

![그림6](https://user-images.githubusercontent.com/49300659/63860495-34f97180-c9e4-11e9-873d-3f5c69b9ea1c.png)


# Requirement

The CoreSampler requires python 3.0 and numpy library.

It also works in python 2.0, but python 3.0 is recommended for handling a large data set.


# Contact

jaeyoonkim72@gmail.com


# Citation

Jeong S, Kim JY, Jeong SC, Kang ST, Moon JK, et al. (2017) GenoCore: A simple and fast algorithm for core subset selection from large genotype datasets. PLOS ONE 12(7): e0181420. https://doi.org/10.1371/journal.pone.0181420
