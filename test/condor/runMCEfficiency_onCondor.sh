#!/bin/bash

cp /afs/cern.ch/user/m/mgaisdor/.globus/x509_proxy x509up_u160839 
export X509_USER_PROXY=x509up_u160839 
export HERE=$PWD
export WORKDIR=${CMSSW_BASE}/src/Analysis/MuonDST-Framework/
export DATASET=$1
export OUTPUT=$2
#export LOGS="/eos/user/${USER:0:1}/$USER/Run3ScoutingMuonReco/logs/"


mkdir -p ${HERE}/logs
#mkdir -p ${LOGS}
mkdir -p ${OUTPUT}

condor_submit runMCEfficiency_onCondor.sub
