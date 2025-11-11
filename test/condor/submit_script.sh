#!/bin/bash

sh runMCEfficiency_onCondor.sh CMSSW_15_1_0/"$1"/DisplacedIterative/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/CMSSW_15_1_0/"$1"/DisplacedIterative/
sh runMCEfficiency_onCondor.sh CMSSW_15_1_0/"$1"/StandardIterative/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/CMSSW_15_1_0/"$1"/StandardIterative/
sh runMCEfficiency_onCondor.sh CMSSW_15_1_0/"$1"/Cascade/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/CMSSW_15_1_0/"$1"/Cascade/