#!/bin/bash

for CLUSTER in IC_2391 Collinder_135 NGC_2451A NGC_2547 IC_2602:
do
    python tess_young/tau_sq/tau_sq_run.py -c config/tau_sq1.yml -g $CLUSTER
    python tess_young/tau_sq/tau_sq_run.py -c config/tau_sq4.yml -g $CLUSTER

    python tess_young/tau_sq/tau_sq_run.py -c config/tau_sq1_log.yml -g $CLUSTER
    python tess_young/tau_sq/tau_sq_run.py -c config/tau_sq4_log.yml -g $CLUSTER
done
