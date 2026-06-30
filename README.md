# FRANZ

FRANZ (FRamework for ANalytical one-Zone blastwave dynamics) is a Julia package for modelling the evolution of astrophysical blastwaves using analytical and semi-analytical methods.

The framework is designed primarily for studying supernova remnants (SNRs) evolving in complex galactic environments. Rather than performing computationally expensive hydrodynamical simulations, FRANZ models blastwave evolution using a one-zone ordinary differential equation (ODE) approach based on the thin-shell and sector approximations. Each angular sector evolves independently, allowing different parts of the blastwave to experience different environmental conditions while remaining computationally inexpensive.

FRANZ aims to bridge the gap between simple analytical solutions and full numerical simulations. It enables rapid parameter studies, facilitates physical interpretation, and provides a modular platform for incorporating additional physical processes.

## Features

* Analytical one-zone model for blastwave evolution
* Thin-shell and sector approximations
* Modular framework for incorporating additional physics
* Rapid exploration of parameter space
* Minimal package dependencies
* Designed for both research and teaching

## Installation

### Installing Julia

If you do not already have Julia installed, download the latest stable release from

https://julialang.org/downloads/

After installation, verify that Julia is available by running

```bash
julia --version
```

For users new to Julia, the official documentation is available at

https://docs.julialang.org/

### Installing FRANZ

**Disclaimer**: The workflow described below is currently not yet supported. Until the package has been registered with the official Julia registries, please simply `git clone` this package to install it.

Once Julia is installed, start a Julia session and install the package:

```julia
using Pkg
Pkg.add("FRANZ")
```

Then load the package with

```julia
using FRANZ
```

## Example

A complete introductory example demonstrating the basic workflow is available in

```text
example/examples.ipynb
```

The notebook introduces the main components of the framework and demonstrates how to construct and evolve a blastwave model.

## Using FRANZ from Python

Although FRANZ is written in Julia, it can be used directly from Python via the `juliacall` package.

First install Julia (see above), then install `juliacall`:

```bash
pip install juliacall
```

A simple Python session might look like

```python
from juliacall import Main as jl

jl.seval("using Pkg")
jl.seval('Pkg.add("FRANZ")')   # only required the first time

jl.seval("using FRANZ")
```

Julia functions and types can then be accessed through `jl`. Since the computational work is still performed in Julia, Python users benefit from the same performance as native Julia users.

## Documentation

The repository currently includes:

* An introductory Jupyter notebook in `example/`
* Inline documentation for the public API

Additional documentation and worked examples will be added as the package develops.  

Basic usage might look like

```Julia
using FRANZ

my_model       = Model(E_51=E_51, M_ej=M_ej, f_CR=f_CR, α_p=α_p, Λ_22=Λ_22)
my_environment = Environment(density=(x::Vector{<:Real}, t::Real) -> crazy_density_field(x, t), 
			     g_ext=(x::Vector{<:Real}, t::Real) -> crazy_gravity_field(x, t))

t, x, v, dA, M, f = numerical_solution(output_times, model=my_model, environment=my_environment, cosθ=0.0, ϕ=LinRange(0, 2π, 256))
```
which defines an explosion model and an environment with a non-standard density and gravitational field and runs the blastwave model for 256 shock-surface segments in the cosθ=0.0 plane.
The output are then time-series of position, velocity, surface normal and mass (per unit solid-angle) for each direction. 

## Contributing

Bug reports, feature requests, and pull requests are welcome. If you encounter unexpected behaviour or would like to suggest new functionality, please open an issue on GitHub.

## Citation

If you use FRANZ in published work, please cite the accompanying publication:

[![DOI](https://zenodo.org/badge/1258322843.svg)](https://doi.org/10.5281/zenodo.21076525)

## License

FRANZ is released under the MIT License.

