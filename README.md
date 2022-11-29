# ARSP

This project require libaraies Qhull and glpk.

## File Description


## Generate Syntehtic Datasets

## Run Algorithms for ARSP

```
./main op path dim m cnt l p c
```
|Parameters|Meaning|Example|
|op|operator|-brachbound, -gendata, ..., (see src-arsp/main.cpp)|
|path|data path|data/inde/, data/anti/|
|dim|data dimensionality|3, 4, ...|
|m|number of uncertain objects|16000, ...|
|cnt|number of instances per object|400, ...|
|l|length of possible region of each object|0.2, ...|
|p|percentage of objects with the empty instance|0.1, ...|
|c|number of constraints|3,4,5,...|