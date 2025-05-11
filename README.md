# PhysicalOptics (DEPRECATED)
A package for simulation of physical optics. Physical optics is more general than ray optics but not as general as full electrodynamics.

> [!WARNING]
> ⚠️⚠️ This package has been deprecated in favor of [WaveOpticsPropagation.jl](https://github.com/JuliaPhysics/WaveOpticsPropagation.jl) ⚠️⚠️


| **Documentation**                       | **Build Status**                | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][CI-img]][CI-url] | [![][codecov-img]][codecov-url] |


## Installation
Currently not registered and under development. But the main branch can be installed with:
```julia
julia> ] add https://github.com/JuliaPhysics/PhysicalOptics.jl
```

## Features
### Implemented
* Fast convolutions adapted to the needs of Physical Optics. Methods like `conv_psf` and `conv_otf` are wrappers for the more general `conv` and `conv_v_ft`. 
* Light propagation with Fresnel and Rayleigh-Sommerfeld (`propagate`) 
* 2D jinc PSF
* some conversion methods
* Light propagation with Fresnel
* Optical elements like lenses, apertures
* some optical conversions
* micro lenses
* Light propagation with Fraunhofer 
* Focused and defocused PSFs (3D) with Debye integral.

### Wanted
* Register adjoints of convolution via ChainRulesCore
* more tests
* cleaning of method arguments
* Documentation

## Literature
As resources we recommend 
* Goodman, Joseph W. Introduction to Fourier optics, 
* Mertz, Jerome. Introduction to optical microscopy. Cambridge University Press, 2019.

For simulation there exists a MATLAB tutorial
* Voelz, David. "Computational fourier optics: a MATLAB tutorial." Society of Photo-Optical Instrumentation Engineers, 2011.


## Similar Projects
In Julia there is no similar project. However, in Python [POPPY](https://github.com/spacetelescope/poppy) offers similar functionality.



[docs-dev-img]: https://img.shields.io/badge/docs-dev-orange.svg 
[docs-dev-url]: https://juliaphysics.github.io/PhysicalOptics.jl/dev/ 

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg 
[docs-stable-url]: https://juliaphysics.github.io/PhysicalOptics.jl/stable/

[codecov-img]: https://codecov.io/gh/JuliaPhysics/PhysicalOptics.jl/branch/main/graph/badge.svg?token=H94RIVDYK4 
[codecov-url]: https://codecov.io/gh/JuliaPhysics/PhysicalOptics.jl 

[CI-img]: https://github.com/JuliaPhysics/PhysicalOptics.jl/workflows/CI/badge.svg
[CI-url]: https://github.com/JuliaPhysics/PhysicalOptics.jl/actions?query=workflow%3ACI 
