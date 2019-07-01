#! /usr/bin/env bash

##################################################

RNDSEED=0
MODE="neutrino"     # Beam mode

GEOMETRY="FullModel"

USERDIR="/dune/app/users/marshalc/ND_neutron"
TOPVOL="volWorld"

##################################################

## Set up UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup dk2nu        v01_05_01b -q e15:prof
setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup ifdhc

##################################################

## Copy some flux files to local folder
./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60

##################################################

## Copy geometry and configuration files to local folder

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

## Run GENIE generating a max path lengths file for the geometry

gevgen_fnal \
    -f local_flux_files/gsimple*.root,DUNEND \
    -g ${GEOMETRY}.gdml \
    -t ${TOPVOL} \
    -m "+${GEOMETRY}.${TOPVOL}.maxpl.xml" \
    -S "+100000" \
    -L cm -D g_cm3 \
    -n 1 \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC

