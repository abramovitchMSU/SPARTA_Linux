#!/bin/bash

#Install Python development headers, NumPy, Java, and R
sudo apt-get update
sudo apt-get install build-essential python2.7-dev python-numpy default-jre r-base r-base-dev
sudo R --vanilla --slave < install_edger.r
