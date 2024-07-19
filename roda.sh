#!/bin/bash

for counter in $(seq 1 10)
do
  /usr/bin/time -a -o timeSerial.txt ./main | tee -a timeSerial.txt

done
