# Sunflow : Efficient Optical Circuit Scheduling for Coflows #

Sunflow is an optical circuit scheduler, with near-packet-switching performance 
for Coflows. Sunflow is described in:

>[Sunflow: Efficient Optical Circuit Scheduling for Coflows](http://www.cs.rice.edu/~eugeneng/papers/CoNEXT16.pdf)

>[Xin Sunny Huang](http://www.cs.rice.edu/~xinh/), 
[Xiaoye Steven Sun](http://www.owlnet.rice.edu/~xs6/), 
and [T. S. Eugene Ng](http://www.cs.rice.edu/~eugeneng/)  

>_In 12th International Conference on emerging Networking EXperiments and 
Technologies ([CoNEXT 2016](http://conferences2.sigcomm.org/co-next/2016/#!/home))_

This repository contains a flow-level, discrete-event simulator for Sunflow, 
along with the tools necessary to replicate the experimental results reported 
in the paper. 

## Zero dependency: run wherever you like ##
The minimal working set is in `src/`, with a `Makefile` to build the simulator 
binary. You can quickly start the simulator with 

``` 
cd src/ 
make 
./Ximulator 
```

Take a look at `main.cc` for the flag options that you may want to reconfigure. 
For example, to play with a specific scheduler, use

``` 
./Ximulator -s sunflow 
```

## IDE ##

This whole package works with CLion, so it comes with lots of `CMakeLists.txt`. 
To use IDE which may have a different directory for the binary other than 
`src/`, you may want to reconfigure the `SRC_DIR` in `global.cc`, 
to help the binary to find out the correct directory for input and output files.
For details, take a look at `global.cc`.

## Tests
I have provided simple tests for all schedulers in 
`tests/test_src/ximulator_test.cc`. To run the test, you may also need to 
reconfigure `TEST_DATA_DIR` in `ximulator_test.cc`. Due to randomization
in the simulation, the tests may fail when running in a different environment.  

## Trace ##
The Coflow trace comes from the popular
[coflow-benchmark](https://github.com/coflow/coflow-benchmark) by 
[Mosharaf Chowdhury](https://github.com/mosharaf), with minor reformatting.   
