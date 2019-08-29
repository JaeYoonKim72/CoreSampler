#!/bin/env python
#Genocore python version 1.0 (made by JaeYoon Kim)
#This is a python version coverted from the orignial R version made by Seongmun Jeong 
#2019-06-01

import os, sys, collections, time, copy
import itertools
import numpy as np
import optparse
from CoreConverter import *
from CoreSampler import *

def parse_my_command_line():
    opt_parser = {}
    opt_parser['VCFtoCSV'] = optparse.OptionParser()
    opt_parser['CoreSampler']  = optparse.OptionParser()
    opt_parser['SelectVCF'] = optparse.OptionParser()

    # VCFtoCSV
    opt_parser['VCFtoCSV'].usage = "python run_CoreSampler.py VCFtoCSV -i [VCF file] -o [Output name] -p [Y or N (phased)] -g [Y or N (gziped)]\n    Help: python genocore.py VCFtoCSV -h"
    opt_parser['VCFtoCSV'].add_option('-i', '--input', default = None, type = 'str', dest='input', action='store',
        help = 'Input vcf file')
    opt_parser['VCFtoCSV'].add_option('-o', '--output', default = "testVCFtoCSV.csv", type = 'str', dest='output', action='store',
        help = 'Output_name.csv')
    opt_parser['VCFtoCSV'].add_option('-p', '--pased', default = 'N', type = 'str', dest='phased', action='store',
        help = 'Wheter the VCf is phased. Select "Y" or "N" (default:"N")')
    opt_parser['VCFtoCSV'].add_option('-g', '--gziped', default = 'N', type = 'str', dest='gziped', action='store',
        help = 'Wheter the VCf is gziped. Select "Y" or "N" (default:"N")')


    # CoreSampler
    opt_parser['CoreSampler'].usage = "python run_CoreSampler.py CoreSampler -i [CSV file] -p [Preset txt] -n [The number of Sample set] -m [MAF rate] -o [Output name] \n    Help: python genocore.py CoreSampler -h"
    opt_parser['CoreSampler'].add_option('-i', '--input', default = None, type = 'str', dest='input', action='store',
        help = 'Input csv file')
    opt_parser['CoreSampler'].add_option('-p', '--preset', default = None, type = 'str', dest='preset', action='store',
        help = 'Preset text file (one sample per line), This is an optional argument and can be omitted if there is no Preset')
    opt_parser['CoreSampler'].add_option('-n', '--number', default = None, type = 'int', dest='nsamples', action='store',
        help = 'The number of samples (default: None)')
    opt_parser['CoreSampler'].add_option('-m', '--maf', default = 0, type = 'float', dest='maf', action='store',
        help = 'Minor allele frequency cut-off (default: 0.0)')
    opt_parser['CoreSampler'].add_option('-o', '--output', default = 'test', type = 'str', dest='output_name', action='store',
        help = 'Output name prefix (default: test)')

    # SelectVCF
    opt_parser['SelectVCF'].usage = "python run_CoreSampler.py SelectVCF -i [VCF file] -g [Y or N (gziped)] -s [Sample list] -o [Output name] \n    Help: python genocore.py VCFtoCSV -h"
    opt_parser['SelectVCF'].add_option('-i', '--input', default = None, type = 'str', dest='input', action='store',
        help = 'Input vcf file')
    opt_parser['SelectVCF'].add_option('-g', '--gziped', default = 'N', type = 'str', dest='gziped', action='store',
        help = 'Wheter the VCf is gziped. Select "Y" or "N" (default:"N")')
    opt_parser['SelectVCF'].add_option('-s', '--sample', default = None, type = 'str', dest='sample', action='store',
        help = 'Sample text file (one sample per line)')
    opt_parser['SelectVCF'].add_option('-o', '--output', default = 'Selected.vcf', type = 'str', dest='output_name', action='store',
        help = 'Output selected vcf file name')


    # Get the command from sys.argv
    try:
        command = sys.argv[1]
    except IndexError:
        command = None

    # Get the appropriate option parser
    try:
        parser = opt_parser[command]
        # Done. Parse arguments and return.
        opts, args = parser.parse_args()
        return parser.parse_args()
    except KeyError:
        # Invalid command. Create a parser to show default usage
        parser       = optparse.OptionParser()
        parser.usage  = '\n\n%prog [Methods (VCFtoCSV or CoreSampler or SelectVCF) ] [options]\n\n'
        parser.usage += 'Description:\n'
        parser.usage += '  This prgram works Python version 3 and requires Numpy library.\n'
        parser.usage += '  To work this program, select Methods, VCFtoCSV, CoreSampler or SelectVCF.\n\n'
        parser.usage += 'Commands for each methods:\n'
        for cmd in ['VCFtoCSV', 'CoreSampler', 'SelectVCF']:
            parser.usage += '  %s\n' % cmd
            parser.usage += '    %s\n\n' % opt_parser[cmd].usage
        parser.usage = parser.usage.strip()
        if command is None:
            parser.error('command cannot be empty')
            sys.exit()
        else:
            parser.error('invalid command')
            sys.exit()

def main():
    options, args = parse_my_command_line()
    if len(args) != 1:
        print(os.system('python run_CoreSampler.py -h'))
        sys.exit()
    return options, args

if __name__ == '__main__':
    opt, args = main()
    if args[0] == 'VCFtoCSV':
        if opt.input == None:
            print(os.system('python run_CoreSampler.py VCFtoCSV -h'))
            sys.exit()
        infile = opt.input
        output = opt.output
        phased = opt.phased
        gziped = opt.gziped
        result = VCF_to_CSV(infile, output, phased, gziped)
        if result == "Done": 
            print("Coverting complete!")
        else:
            print("Coverting Error")
            sys.exit()

    if args[0] == "CoreSampler":
        if opt.input == None or opt.nsamples == None:
            print(os.system('python run_CoreSampler.py CoreSampler -h'))
            sys.exit()
        infile = opt.input
        preset = opt.preset
        outname = opt.output_name
        maf = opt.maf 
        nsample = opt.nsamples
        result = CoreSampler(infile, preset, nsample, maf, outname)
        if result == "Done":
            print("CoreSampler completed!")
        else:
            print("CoreSampler incompleted!")
            sys.exit()
  
    if args[0] == 'SelectVCF':
        if opt.input == None or opt.sample == None:
            print(os.system('python run_CoreSampler.py SelectVCF -h'))
            sys.exit()
        infile = opt.input
        gziped = opt.gziped
        sample = opt.sample
        outname = opt.output_name
        result = VCF_select(infile, sample, outname, gziped)
        if result == "Done":
            print("VCF selection completed!")
        else:
            print("VCF selection incompleted!")
            sys.exit()
