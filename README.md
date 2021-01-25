# Schistoxpkg

Go to https://mattg3004.github.io/Schistoxpkg.jl/dev/ for documentation

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add Schistoxpkg
```
Alternatively directly from Julia REPL run:
```
julia> using Pkg
julia> Pkg.add("Schistoxpkg")
```

Additionally, there is an issue with the save and load package JLD, which will produce an error if you have the latest versions of all libraries installed. To get around this problem, you need to install a previous version HFD5_jll library as follows:
```
julia> Pkg.add(Pkg.PackageSpec(;name="HDF5_jll", version=“1.10.5”))
```

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mattg3004.github.io/Schistoxpkg.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mattg3004.github.io/Schistoxpkg.jl/dev)

[![Build Status](https://travis-ci.com/mattg3004/Schistoxpkg.jl.svg?branch=master)](https://travis-ci.com/mattg3004/Schistoxpkg.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mattg3004/Schistoxpkg.jl?svg=true)](https://ci.appveyor.com/project/mattg3004/Schistoxpkg-jl)

[![Codecov](https://codecov.io/gh/mattg3004/Schistoxpkg.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mattg3004/Schistoxpkg.jl)


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg



