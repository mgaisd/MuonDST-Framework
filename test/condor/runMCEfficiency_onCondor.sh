#!/bin/bash

cp /tmp/x509up_u160839 x509up_u160839 
export X509_USER_PROXY=x509up_u160839 
export HERE=$PWD
export WORKDIR=${CMSSW_BASE}/src/Analysis/MuonDST-Framework/
export DATASET=$1
export OUTPUT=$2

mkdir -p ${HERE}/logs/
mkdir -p ${OUTPUT}

condor_submit runMCEfficiency_onCondor.sub
