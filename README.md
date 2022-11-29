# ARSP

This project require libaraies Qhull and glpk.

## File Description


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

### generate F

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
./main -baseline-V/-branchbound/-kdtree/-kdtree-star/-quadtree-star path dim m cnt l p c
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


### run src-eclipse