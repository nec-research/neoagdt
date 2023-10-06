# neoag-digital-twin
A digital twin software framework for optimal neoantigen-based treatments.

### PIP installation
`pip install .[all]`

### Conda installation
Alternatively, we provide an Anaconda-based installation:

```
conda env create -f environment.yml  # create environment with all dependencies
conda activate neoagdt
pip install .  # install the neoagdt Python package (dependencies should already be insalled at the previous step)
```

Python 3.9 is required.

---

### Sample configurations
Configurations are in `etc/`.
Paths can also be absolute and do not need to be relative to the repository.

### Simulate a cell population
```
simulate-cancer-cells etc/cells-config.yaml --logging-level INFO
```

### Run vaccine optimization
The vaccine design optimization supports multiprocessing leveraging on a local `dask` cluster.
It distributes calls to the MIP optimization. 
* `--num-procs <NUM_PROCS>`: number of parallel processes (and of CPU cores)
* `--num-threads-per-proc <NUM_THREADS_PER_PROC>`: number of threads to allocate for each process. So the total number of threads for a local cluster will be (`args.num_procs * args.num_threads_per_proc`)
```
optimize-vaccine-ilp etc/optimization-config.yaml --num-procs <NUM_PROCS> --num-threads-per-proc <NUM_THREADS_PER_PROC> --logging-level INFO
```

### Plot bar charts
```
create-bar-chart etc/bar-charts-config.yaml --logging-level INFO
```

### Evaluate vaccine response
```
evaluate-vaccine-response etc/response-likelihood-config.yaml --logging-level INFO
```
---

### Execute unit tests
```
pip install .[test]
pytest tests
```

### Execute style tests
```
pip install .[test]
pylint neoag_dt
```

### Documentation

The documentation project for this project can be built with `sphinx`. The necessary dependencies are install by pip
when installing either the `all` or `docs` optional dependencies.

```
pip install .[docs]
cd docs
make html
```