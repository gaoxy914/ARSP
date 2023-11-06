# ARSP

This project require c++ libraries Qhull and glpk and python libraries nba_api, pandas, and numpy.

## File Description

|File|Description|
|:---|:---|
|IIP.py|Generating IIP from iip_98_2000.txt|
|CAR.py|Generating CAR from autos.csv|
|NBA.py|For getting nba dataset from https://www.nba.com/stats/.|
|src-arsp|source code of algorithms for ARSP. All avaiable methods can be found in all_rsky_prob.h and region.h|
|src-eclp|source code of algorihtms for ECLP. All avaiable methods can be found in all_ecl_prob.h|
|src-eclipse|source code of algorithms for ECLIPSE. All avaiable methods can be found in eclipse.h|

## Commands for Complation and Execution

### complie

Change SRC_DIR in Makefile to the project you want to complie and run `make`.

### generate nba dataset

Run `python nba_dataset.py` to get file nba.data and nba.cnt.

### generate synthetic datasets

After compiling the src-arsp, run
```
./main -gendata path dim m cnt l p c
```
|Parameters|Meaning|Example|
|:---:|:---|:---:|
|path|data path|data/inde/, data/anti/|
|dim|data dimensionality|3, 4, ...|
|m|number of uncertain objects|16000, ...|
|cnt|number of instances per object|400, ...|
|l|length of possible region of each object|0.2, ...|
|p|percentage of objects with the empty instance|0.1, ...|
|c|number of constraints|3,4,5,...|

### generate the scoring function set F

After compiling the src-arsp, run
```
./main -genquery dim c
```
|Parameters|Meaning|Example|
|:---:|:---|:---:|
|dim|data dimensionality|3, 4, ...|
|c|number of constraints|3,4,5,...|

### run src-arsp

After compiling the src-arsp, run
```
./main -enum/-baseline-V/-branchbound/-kdtree/-kdtree-star/-quadtree-star path dim m cnt l p c
```
|Parameters|Meaning|Example|
|:---:|---|:---:|
|path|data path|data/inde/, data/anti/|
|dim|data dimensionality|3, 4, ...|
|m|number of uncertain objects|16000, ...|
|cnt|number of instances per object|400, ...|
|l|length of possible region of each object|0.2, ...|
|p|percentage of objects with the empty instance|0.1, ...|
|c|number of constraints|3,4,5,...|

### run src-eclp

After compiling the src-eclp, run
```
./main -dua-ms/-kdtree-star/ path dim m cnt l p
```
|Parameters|Meaning|Example|
|:---:|---|:---:|
|path|data path|data/inde/, data/anti/|
|dim|data dimensionality|3, 4, ...|
|m|number of uncertain objects|16000, ...|
|cnt|number of instances per object|400, ...|
|l|length of possible region of each object|0.2, ...|
|p|percentage of objects with the empty instance|0.1, ...|

### run src-eclipse

After compiling the src-eclipse, run `./main dim n` to compare the performance of DUAL-S and QUAD (the state-of-the-art for eclipse query processing).

|Parameters|Meaning|Example|
|:---:|---|:---:|
|dim|data dimensionality|3, 4, ...|
|n|number of tuples|10000, ...|