# MomentArithmetic.jl

Rigorous propagation of partial inforation about [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) and [dependencies](https://en.wikipedia.org/wiki/Copula_(probability_theory)) in Julia.

Generalisation of First Order Error Propagation with:

- Independence between variables need not be assumed (but can be)
- Moments may be partially known (intervals)
- Input distribution assumptions no longer necessary
- Distributrional information can be obtained from moment and range information in non-linear models

May be viewed as a form of *distribution-free risk analysis*.

[REC21 Paper](https://elements.liverpool.ac.uk/repository.html?com=get-file&publication-id=495120&rfid=https%3A%2F%2Flivrepository.liverpool.ac.uk%2Frt4eprints%2Ffile%2F393269%2FREC2021_37_Gray.pdf)

[REC21 Presentation](http://ww2new.unime.it/REC2021/index.php?uri=presentations)


## Installation

Not yet a registered package, however may be installed through the Julia package manager:

```Julia
julia> ]
pkg> add https://github.com/AnderGray/MomentArithmetic.jl
```


## Features

- 9 univariate transformations: scalar multiplication, scalar translation, exp, log, ln, 1/x, x^2, sqrt and |x|
- 7 binary operations: +, -, \*, /, ^, min, max
- Independence and no dependence assumptions (Fréchet)
- Moments consitancy checking
- Bounding [p-box](https://en.wikipedia.org/wiki/Probability_box) from Moments

## Usage

TODO
