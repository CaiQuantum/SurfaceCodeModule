# SurfaceCodeModule
This is a python module written in C++ that can perform threshold simulation, using Blossom V minimum weight perfect matching as decoder.
## Build the module
Download the repo which will be in a folder called `SurfaceCodeModule-master`. Open your terminal, change your directory to the `SurfaceCodeModule-master` folder and input
```bash
python setup.py build
```
A .so python module file will be built in the `SurfaceCodeModule-master/build` folder, which can be imported into python.

## Usage
Assuming the .so file is located in module_path, we can use it in the following way:

```python
import sys
sys.path.append(module_path)
from SurfaceCode import average_logical_error_array
log_err_array = average_logical_error_array(code_distance, n_runs, X_error_table, Z_error_table, 
                                            planar_switch, time_step_ratio, time_weight_ratio)
```
Here `average_logical_error_array` is the function we use to calculate the average logical error rate of the code given the following arguments:
* `code_distance`: distance of the surface code.
* `n_runs`: Number of sample runs
* `X_error_table`, `Z_error_table`: the error tables that outline the proability of different Pauli and measurement error probabilities for X stabiliser checks and Z stabiliser checks circuits.
* `planar_switch`: 0 means toric code, 1 means planar code.
* `time_step_ratio`: the number of rounds of repeated parity measurements versus the code distance. 
* `time_weight_ratio`: the weight of a edge in time direction versus that in spatial direction when performing minimum weight perfect matching. 

The returned list `log_err_array` consist of 5 elements: 
```python
[total_logical_error_rate, X1_logical_error_rate, Z1_logical_error_rate, X2_logical_error_rate, Z2_logical_error_rate]
```
here 1 and 2 denote the label of the logical qubits. For toric code, we have two logical qubits, thus the whole array is filled. For planar code however, we have only one logical qubit, thus the last two arguments the list are meaningless when we are considering planar code.


