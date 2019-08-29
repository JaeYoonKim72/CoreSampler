# CoreSampler: core samples/subset selection from large genotype datasets

CoreSampler is a modified version of GenoCore, which improves on the limitations of GenoCore where the user is unable to specify the number of samples to extract.

The CoreSampler extracts core samples as many numbers as user wanted, and provides a VCF file of core samples and a coverage file recorded how well the core samples are representative to the entire sample.

This program is based on the GenoCore python version availiable at https://github.com/JaeYoonKim72/GenoCore_Python.

Source code was written by Jae-Yoon Kim using Python language and supported on windows and linux platform.


# CoreSampler algorithm

The "CoreSampler" consists of three submodules: VCFtoCSV, CoreSampler, SelectVCF

"CoreSampler", the main submodule, sets the number of samples to extract via the "-n" argument.

The algorithm schematic is as follows.

![그림6](https://user-images.githubusercontent.com/49300659/63860495-34f97180-c9e4-11e9-873d-3f5c69b9ea1c.png)


# Usage

## 1) Download

git clone https://github.com/JaeYoonKim72/CoreSampler

cd GenoCore_Python


## 2) Basic usage

Basic usage : run_CoreSampler.py [Methods (VCFtoCSV or CoreSampler or SelectVCF) ] [options]

![basic](https://user-images.githubusercontent.com/49300659/63917240-25763900-ca75-11e9-8bbf-2603ae7d0109.png)


## 3) VCFtoCSV
The "VCFtoCSV" submodule converts a VCF file into a CSV file for CoreSampler. The output file is a CSV file for that VCF.

Usage: python run_CoreSampler.py VCFtoCSV -i [VCF file] -o [Output name] -p [Y or N (phased)] -g [Y or N (gziped)]

Example: python run_CoreSampler.py VCFtoCSV -i ExampleData/Test_420sample.vcf.gz -o ExampleData/Test_420sample.csv -p Y -g Y

![VCFtoCSV](https://user-images.githubusercontent.com/49300659/63917251-2e670a80-ca75-11e9-8bd6-9d3dd93e5222.png)

## 4) CoreSampler
The "CoreSampler" submodule is the main module and performs core sample extraction. 

The number of samples to be extracted is set via the "-n" argument, and the "Preset" file, specifying samples that must be extracted, is set via the "-p" argument. This argument can be omitted if there are no preset samples.

The output files are core-sample-list, core-sample-csv and core-sample-coverage files, and removed-marker file due to MAF.

Usage: python run_CoreSampler.py CoreSampler -i [CSV file] -p [Preset txt] -n [The number of Sample set] -m [MAF rate] -o [Output name]

Example: python run_CoreSampler.py CoreSampler -i ExampleData/Test_420sample.csv -p ExampleData/Preset.txt -n 30 -m 0.05 -o ExampleData/TestCoreSet

![CoreSampler](https://user-images.githubusercontent.com/49300659/63917265-3626af00-ca75-11e9-9645-59a2395a33ba.png)

## 5) SelectVCF
The "SelectVCF" submodule extract a core-set VCF file, using the core-sample-list file created by the "CoreSampler" module and the input VCF file used by "VCFtoCSV". The output file is a VCF file for final core samples.

Usage: python run_CoreSampler.py SelectVCF -i [VCF file] -g [Y or N (gziped)] -s [Sample list] -o [Output name]

Example: python run_CoreSampler.py SelectVCF -i ExampleData/Test_420sample.vcf.gz -g Y -s ExampleData/TestCoreSet_CoreSample.list.txt -o ExampleData/TestCoreSet_CoreSampler.vcf

![SelectVCF](https://user-images.githubusercontent.com/49300659/63917279-3de65380-ca75-11e9-9d6c-679830e4676b.png)


# Requirement

The CoreSampler requires python 3.0 and numpy library.

It also works in python 2.0, but python 3.0 is recommended for handling a large data set.


# Contact

jaeyoonkim72@gmail.com


# License

CoreSampler is registered with the Korean Copyright Commission under accession umber C-2017-024343.


# Citation

Jeong S, Kim JY, Jeong SC, Kang ST, Moon JK, et al. (2017) GenoCore: A simple and fast algorithm for core subset selection from large genotype datasets. PLOS ONE 12(7): e0181420. https://doi.org/10.1371/journal.pone.0181420
