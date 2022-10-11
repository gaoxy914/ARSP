#!/bin/bash
echo "m = 4096" >> time-inde.txt
./trans data/inde/ 4096 >> time-inde.txt
./trans data/inde/ 4096 >> time-inde.txt
./trans data/inde/ 4096 >> time-inde.txt
echo "m = 16384" >> time-inde.txt
./trans data/inde/ 16384 >> time-inde.txt
./trans data/inde/ 16384 >> time-inde.txt
./trans data/inde/ 16384 >> time-inde.txt
echo "m = 65536" >> time-inde.txt
./trans data/inde/ 65536 >> time-inde.txt
./trans data/inde/ 65536 >> time-inde.txt
./trans data/inde/ 65536 >> time-inde.txt
echo "m = 262144" >> time-inde.txt
./trans data/inde/ 262144 >> time-inde.txt
./trans data/inde/ 262144 >> time-inde.txt
./trans data/inde/ 262144 >> time-inde.txt

echo "m = 256" >> time-anti.txt
./trans data/anti/ 256 >> time-anti.txt
./trans data/anti/ 256 >> time-anti.txt
./trans data/anti/ 256 >> time-anti.txt
echo "m = 1024" >> time-anti.txt
./trans data/anti/ 1024 >> time-anti.txt
./trans data/anti/ 1024 >> time-anti.txt
./trans data/anti/ 1024 >> time-anti.txt
echo "m = 4096" >> time-anti.txt
./trans data/anti/ 4096 >> time-anti.txt
./trans data/anti/ 4096 >> time-anti.txt
./trans data/anti/ 4096 >> time-anti.txt
echo "m = 16384" >> time-anti.txt
./trans data/anti/ 16384 >> time-anti.txt
./trans data/anti/ 16384 >> time-anti.txt
./trans data/anti/ 16384 >> time-anti.txt
echo "m = 65536" >> time-anti.txt
./trans data/anti/ 65536 >> time-anti.txt
./trans data/anti/ 65536 >> time-anti.txt
./trans data/anti/ 65536 >> time-anti.txt
echo "m = 262144" >> time-anti.txt
./trans data/anti/ 262144 >> time-anti.txt
./trans data/anti/ 262144 >> time-anti.txt
./trans data/anti/ 262144 >> time-anti.txt

echo "m = 256" >> time-corr.txt
./trans data/corr/ 256 >> time-corr.txt
./trans data/corr/ 256 >> time-corr.txt
./trans data/corr/ 256 >> time-corr.txt
echo "m = 1024" >> time-corr.txt
./trans data/corr/ 1024 >> time-corr.txt
./trans data/corr/ 1024 >> time-corr.txt
./trans data/corr/ 1024 >> time-corr.txt
echo "m = 4096" >> time-corr.txt
./trans data/corr/ 4096 >> time-corr.txt
./trans data/corr/ 4096 >> time-corr.txt
./trans data/corr/ 4096 >> time-corr.txt
echo "m = 16384" >> time-corr.txt
./trans data/corr/ 16384 >> time-corr.txt
./trans data/corr/ 16384 >> time-corr.txt
./trans data/corr/ 16384 >> time-corr.txt
echo "m = 65536" >> time-corr.txt
./trans data/corr/ 65536 >> time-corr.txt
./trans data/corr/ 65536 >> time-corr.txt
./trans data/corr/ 65536 >> time-corr.txt
echo "m = 262144" >> time-corr.txt
./trans data/corr/ 262144 >> time-corr.txt
./trans data/corr/ 262144 >> time-corr.txt
./trans data/corr/ 262144 >> time-corr.txt
