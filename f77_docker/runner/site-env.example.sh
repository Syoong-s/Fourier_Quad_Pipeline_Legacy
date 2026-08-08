#!/usr/bin/env bash

# This optional file is sourced before HPC_MODULES from f77pipeline.env.
# Use it only when the cluster needs extra shell initialization before the
# standard module command becomes available.
#
# source /etc/profile.d/modules.sh
#
# Configure normal compiler, MPI, and Apptainer modules in HPC_MODULES so both
# the pipeline and MPI smoke jobs use the same environment.
