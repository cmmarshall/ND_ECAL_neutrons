#! /usr/bin/env bash

####################

HORN=$1
NPER=$2
FIRST=$3
TOPVOL=$4
GENIEGEOM=$5
TEST=$6
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

if [ "${GENIEGEOM}" = "" ]; then
GENIEGEOM="FullModel"
fi
echo "Running GENIE with geomtry: ${GENIEGEOM}"

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

# use narrower flux window unless we are doing rock
#if [ "${TOPVOL}" != "World" ]; then
#cp GNuMIFlux_narrow.xml GNuMIFlux.xml
#fi

# Get flux files to local node
chmod +x copy_dune_ndtf_flux
#./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/dbrailsf/flux/nd/gsimple/v2_8_6d --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=60
./copy_dune_ndtf_flux --top /pnfs/dune/persistent/users/ljf26/fluxfiles/g4lbne/v3r5p4/QGSP_BERT --dk2nu --outpath local_flux_files --flavor ${MODE} --base OptimizedEngineeredNov2017 --maxmb=300

####################

## Add the location of the GNuMIFlux.xml to the GENIE xml path

export GXMLPATH=${PWD}:${GXMLPATH}
export GNUMIXML="GNuMIFlux.xml"

## Run GENIE and copy output file to dCache persistent

echo "maxpl file: ${GENIEGEOM}.${VOLVOL}.maxpl.xml"

gevgen_fnal \
    -f local_flux_files/*dk2nu*.root,DUNEND \
    -g ${GENIEGEOM}.gdml \
    -t ${VOLVOL} \
    -m ${GENIEGEOM}.${VOLVOL}.maxpl.xml \
    -L cm -D g_cm3 \
    ${NEVENTS} \
    --seed ${RNDSEED} \
    -r ${RNDSEED} \
    -o ${MODE} \
    --message-thresholds Messenger_production.xml \
    --cross-sections ${GENIEXSECPATH}/gxspl-FNALsmall.xml \
    --event-record-print-level 0 \
    --event-generator-list Default+CCMEC

####################
# edep-sim
####################

cp ${MODE}.${RNDSEED}.ghep.root input_file.ghep.root

setup genie        v2_12_10c  -q e15:prof
setup genie_xsec   v2_12_10   -q DefaultPlusValenciaMEC
setup genie_phyopt v2_12_10   -q dkcharmtau
setup geant4       v4_10_3_p01b -q e15:prof

${CP} ${USERDIR}/edep-sim.tar.gz ${PWD}/edep-sim.tar.gz
tar -xzf edep-sim.tar.gz

gntpc -i input_file.ghep.root -f rootracker \
      --event-record-print-level 0 \
      --message-thresholds Messenger_production.xml

# fix the vertex times of the gtracker file
root -l -q -b hack.C

rm -f input_file.ghep.root

export LD_LIBRARY_PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/lib:${LD_LIBRARY_PATH}
export PATH=${PWD}/edep-sim/edep-gcc-6.4.0-x86_64-pc-linux-gnu/bin:${PATH}

# run all the events in a file
#if [ "${NPER}" = -1 ]; then
#echo "Specified all events, determining how many events there are"
#NPER=$(echo "std::cout << gtree->GetEntries() << std::endl;" | genie -l -b input_file.ghep.root 2>/dev/null  | tail -1)
#echo "There are ${NPER} events"
#fi

##################################################

## Run edep-sim
echo "edep-sim -C -g ${GEOMETRY}.gdml -o ${PWD}/output_file.root -u -e ${NPER} dune-nd.mac"
edep-sim \
    -C \
    -g ${GEOMETRY}.gdml \
    -o ${PWD}/output_file.root \
    -u \
    -e ${NPER} \
    dune-nd.mac

## Copy output file to dCache persistent
echo "${CP} output_file.root ${USERDIR}/EDep/${HORN}/${TOPVOL}/${MODE}.${RNDSEED}.edepsim.root"
${CP} output_file.root ${USERDIR}/EDep/${HORN}/${TOPVOL}/${MODE}.${RNDSEED}.edepsim.root
rm -f output_file.root

echo "${CP} ${MODE}.${RNDSEED}.ghep.root ${USERDIR}/GENIE/${HORN}/${TOPVOL}/${MODE}.${RNDSEED}.ghep.root"
${CP} ${MODE}.${RNDSEED}.ghep.root ${USERDIR}/GENIE/${HORN}/${TOPVOL}/${MODE}.${RNDSEED}.ghep.root
rm -f ${MODE}.${RNDSEED}.ghep.root


