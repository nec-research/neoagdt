# neoag-digital-twin
A digital twin software framework for optimal neoantigen-based treatments.

### Installation
`pip install .[all]`

Requires Python 3.8

---

### Sample configurations
Configurations are in `etc/`.
Paths can also be made absolute and do not need to be relative to the repository.

### Simulate a cell population
```
cd neoag-digital-twin
simulate-cancer-cells etc/cells-config.yaml --logging-level INFO
```

### Run vaccine optimization
The vaccine design optimization supports multiprocessing leveraging on a local `dask` cluster.
It distributes calls to the MIP optimization. 
* `--num-procs <NUM_PROCS>`: number of parallel processes (and of CPU cores)
* `--num-threads-per-proc <NUM_THREADS_PER_PROC>`: number of threads to allocate for each process. So the total number of threads for a local cluster will be (`args.num_procs * args.num_threads_per_proc`)
```
cd neoag-digital-twin
optimize-vaccine-ilp etc/optimization-config.yaml --num-procs <NUM_PROCS> --num-threads-per-proc <NUM_THREADS_PER_PROC> --logging-level INFO
```

### Plot bar charts
```
cd neoag-digital-twin
create-bar-chart etc/bar-charts-config.yaml --logging-level INFO
```

### Evaluate vaccine response
```
cd neoag-digital-twin
evaluate-vaccine-response etc/response-likelihood-config.yaml --logging-level INFO
```
---

### Execute unit tests
```
cd neoag-digital-twin
pytest tests
```

### Execute style tests
```
cd neoag-digital-twin
pylint neoag_dt
```

### Docker
To manually build a Docker image:
```
cd neoag-digital-twin
docker build --tag neoag-digital-twin:<version> .
```

### Documentation

The documentation project for this project can be built with `sphinx`. The necessary dependencies are install by pip
when installing either the `all` or `docs` optional dependencies.

```
cd docs
make html
```