# Stochastic Multistage Optimization - SMPS utilities for Python

Parse datasets stored in the SMPS file format. 

Particular emphasis has been put in maintaining sparse representation of the data, to allow
processing of large instances. Scenario matrices are generated lazily.


## Introduction

[SMPS](http://myweb.dal.ca/gassmann/smps2.htm) is a file format for distributing **Stochastic Multistage Optimization Models**. 
It is an extension of the [MPS file format](https://en.wikipedia.org/wiki/MPS_(format)), which is used to store 
optimization models. 


There is a [parser written in FORTRAN](http://myweb.dal.ca/gassmann/inputs.htm), but I was unable to make it work in my 
Python environment. [Pyomo](http://www.pyomo.org/) has some code to write SMPS files, but not to parse them (see also 
[my request](https://groups.google.com/forum/#!searchin/pyomo-forum/smps/pyomo-forum/jfRD7BK4Mt4/GQqbpAzaBAAJ) on their forum). 
I decided to write a parser script from scratch. 

#### Datasets archives:

This library can be used to parse most models from the following archives:
 * [SIPLIB](http://www2.isye.gatech.edu/~sahmed/siplib/)
 * [POST](http://users.iems.northwestern.edu/~jrbirge/html/dholmes/post.html)
 * [Andy Felt's collection](http://www4.uwsp.edu/math/afelt/slptestset/download.html).

## Installation

~~~~
>> git clone https://github.com/robin-vjc/smps
>> cd smps/
>> pip install -e .
~~~~

Optionally, run tests suite
~~~~
>> py.test
~~~~

## Example usage

* Parse stochastic model from smps dataset:
~~~~
from smps.read import StochasticModel
sm = StochasticModel('./data/invest')
~~~~

Expects the files './data/invest.cor', './data/invest.tim' and './data/invest.sto' to exist and follow the
SMPS format.

* Plot scenario tree:

~~~~
sm.plot_scenario_tree()
~~~~

* Generate deterministic equivalent as a single Gurobi model.

~~~~
det_eq = sm.generate_deterministic_equivalent()
det_eq.optimize()
~~~~

* Attributes in StochasticModel of particular interest


    - StochasticModel.A: 
    
    A dictionary of sparse matrices; each entry of A is a sub-block of the nominal model matrix A. 
    E.g., A[3,2] corresponds to the sub-block of A corresponding to the variables of the second
    period as they enter in the constraints of the third period of the problem.
        
    - b, c, lb, ub, vtype:  
    
    Also dictionaries, one entry in the dictionary per stage, and each entry a vector with the
    nominal problem data.
                        
    - StochasticModel.scenario['SCENARIO_LABEL'].A: 
    
    modified matrix A according to scenario 'SCENARIO_LABEL'. Same way b, c, lb, ub and vtype can be accessed. 
    These are generated lazily.
                                                
    - StochasticModel.generate_deterministic_equivalent():  
    
    Returns a gurobi model for the deterministic equivalent of the problem. Only available for models of type 
    scenarios discrete.

## Supported Model Types

The parser does support `SCENARIO DISCRETE` types (which are by far the most commonly found in datasets archives); 
it does not currently support `BLOCK-INDEP` mode.

Also, while standards are set in the SMPS documentation, certain models archives do not always follow 
them, which may cause parsing problems. Please report them as bugs in this case.
