#!/bin/bash

sh runMCEfficiency_onCondor.sh "$1"/DisplacedIterative/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/v2/"$1"/DisplacedIterative/
sh runMCEfficiency_onCondor.sh "$1"/StandardIterative/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/v2/"$1"/StandardIterative/
sh runMCEfficiency_onCondor.sh "$1"/Cascade/ /eos/user/m/mgaisdor/Run3ScoutingMuonReco/v2/"$1"/Cascade/
