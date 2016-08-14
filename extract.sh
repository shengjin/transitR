#!/bin/bash

#line,Rearth,Pressure,Temp,Mearth,tau,Tint,Teq

awk '{print NR"  "$39"  "$48"  "$49"  "$37"  "0.6667"  "$55"  "$50}' ref_red5e9.dat > LineRearthPTMearthTauTintTeq.dat
