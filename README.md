# Welcome to TomoForward.jl

TomoForward is a Julia wrapper of ASTRA toolbox [1] for tomographic forward projection operators to generate sinogram data from 2D images. (3D will be supported soon) 

ProjGeom data structure constructs a projection geometry and its interface is consistent with ASTRA-toolbox.

## Install

Install [Julia](https://julialang.org/downloads/) and in [Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/),

```
] add https://github.com/JuliaTomo/TomoForward.jl
```

## Examples

Please see codes in `examples` folder.


# Features

## Forward projection

### Forward operator for 2D image as a sparse matrix

- 2D paralleal beam with line and strip projection model. Ported from [1]

# Reference

- [1] https://astra-toolbox.com
