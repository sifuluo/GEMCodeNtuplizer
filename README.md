# GEMCodeNtuplizer

## Purpose

This tool is an addon for GEMCode (https://github.com/gem-sw/GEMCode)
It is to be placed inside of GEMCode and used to ntuplize the output of digitized and matched data samples.

## Instruction

CMSSW and GEMCode is required before setting up this package.
After GEMCode is set up, run the following recipe to set up this package.

    cd $CMSSW_BASE/src/GEMCode/GEMValidation/test
    git init
    git remote add origin git@github.com:sifuluo/GEMCodeNtuplizer.git
    git pull origin master
    scram b -j 4

By default, the configuration file in `Project/` are used for testing local samples.
And you will need a `out/` folder to store the test outputs.

To run batch job, use the `MakeSubmission.py` to create required folders and scripts.
