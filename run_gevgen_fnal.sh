#! /usr/bin/env bash

####################

HORN=$1
NPER=$2
FIRST=$3
TOPVOL=$4
TEST=$5
if [ "${HORN}" != "FHC" ] && [ "${HORN}" != "RHC" ]; then
echo "Invalid beam mode ${HORN}"
echo "Must be FHC or RHC"
kill -INT $$
fi

if [ "${NPER}" = "" ]; then
echo "Number of events per job not specified, using 1000"
NPER=1000
fi

if [ "${FIRST}" = "" ]; then
echo "First run number not specified, using 0"
FIRST=0
fi

if [ "${TOPVOL}" = "" ]; then
echo "Top volume not specified, using World"
TOPVOL="World"
fi
VOLVOL="vol${TOPVOL}"

CP="ifdh cp"
if [ "${TEST}" = "test" ]; then
echo "In TEST mode, assuming interactive running"
PROCESS=0
mkdir -p test
cd test
fi

MODE="neutrino"
if [ "${HORN}" = "RHC" ]; then
MODE="antineutrino"
fi

echo "Running gevgen for ${NPER} events in ${HORN} mode"

RNDSEED=$((${PROCESS}+${FIRST}))
NEVENTS="-n ${NPER}"      # No. of events, -e XE16 for POT

GEOMETRY="FullModel"
USERDIR="/pnfs/dune/persistent/users/marshalc/neutronSim"

####################

## Setup UPS and required products

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup dk2nu        v01_05_01b -q e15:prof
setup ifdhc

####################

## Copy stuff to the local node
${CP} ${USERDIR}/stuff.tar.gz stuff.tar.gz
tar -xzf stuff.tar.gz

mv stuff/* .

# Get flux files to local node
chmod +x copy_dune_ndtf_flux
#./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60
./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT --dk2nu --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60

####################

## Add the location of the GNuMIFlux.xml to the GENIE xml path

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

## Run GENIE and copy output file to dCache persistent

gevgen_fnal \
    -f local_flux_files/*dk2nu*.root,DUNEND \
    -g ${GEOMETRY}.gdml \
    -t ${VOLVOL} \
    -m ${GEOMETRY}.${VOLVOL}.maxpl.xml \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC

if [ ${TEST} != "test" ]; then
        ${CP} ${MODE}.${RNDSEED}.ghep.root ${USERDIR}/GENIE/${HORN}/${TOPVOL}/${MODE}.${RNDSEED}.ghep.root
        rm -f ${MODE}.${RNDSEED}.ghep.root
fi
