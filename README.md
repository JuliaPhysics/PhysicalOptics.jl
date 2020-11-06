# PhysicalOptics
A package for simulation of physical optics. Physical optics is more general than ray optics but not as general as full electrodynamics.

| **Documentation**                       | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-dev-img]][docs-dev-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |


# Installation
Currently not registered and under development. But the current main branch can be installed with:
```julia
julia> ] add https://github.com/JuliaPhysics/PhysicalOptics.jl
```


# Features

## Implemented Features
* Fast convolutions adapted to the needs of Physical Optics. Methods like `conv_psf` and `conv_otf` are wrappers for the more general `conv` and `conv_v_ft`. 

## Wanted Features
* Light propagation with Fresnel, Fraunhofer and Rayleigh-Sommerfeld 
* Optical elements like lenses, apertures, micro lenses
* Focused and Defocused PSFs (2D and 3D). For example with Debye integral.

# Similar Projects
In Julia there is no similar project. However, in Python [POPPY](https://github.com/spacetelescope/poppy) offers similar functionality.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-orange.svg 
[docs-dev-url]: https://juliaphysics.github.io/PhysicalOptics.jl/dev/ 

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg 
[docs-stable-url]: https://juliaphysics.github.io/PhysicalOptics.jl/stable/

[travis-img]: https://api.travis-ci.com/JuliaPhysics/PhysicalOptics.jl.svg?branch=main&status=created 
[travis-url]: https://travis-ci.com/github/JuliaPhysics/PhysicalOptics.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/abxnasacbo42jqvc?svg=true 
[appvey-url]: https://ci.appveyor.com/project/roflmaostc/physicaloptics-jl 

[coveral-img]: https://coveralls.io/repos/github/JuliaPhysics/PhysicalOptics.jl/badge.svg?branch=main
[coveral-url]: https://coveralls.io/github/JuliaPhysics/PhysicalOptics.jl?branch=main

[codecov-img]: https://codecov.io/gh/JuliaPhysics/PhysicalOptics.jl/branch/main/graph/badge.svg?token=H94RIVDYK4 
[codecov-url]: https://codecov.io/gh/JuliaPhysics/PhysicalOptics.jl 
