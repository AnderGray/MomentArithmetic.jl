# MomentArithmetic.jl

Rigorous moment propagation with partial information about [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) and [dependencies](https://en.wikipedia.org/wiki/Copula_(probability_theory)) in Julia.

Generalisation of First Order Error Propagation with:

- Independence between variables need not be assumed (but can be)
- Moments may be partially known (intervals)
- Input distribution assumptions no longer necessary
- Distributional information can be obtained from moment and range information in non-linear models

May be viewed as a form of *distribution-free risk analysis*.

[REC21 Conference Paper](https://www.researchgate.net/publication/352225779_Distribution-free_uncertainty_propagation)

[REC21 Presentation](http://ww2new.unime.it/REC2021/presentations/REC2021-37_Presentation.pdf)

## Features

- 9 univariate transformations: scalar multiplication, scalar translation, exp, log, ln, 1/x, x^2, sqrt and |x|
- 7 binary operations: +, -, \*, /, ^, min, max
- Independence and no dependence assumptions (Fréchet)
- Moments consitancy checking
- Bounding [p-box](https://en.wikipedia.org/wiki/Probability_box) from Moments

## Installation

This is a registered Julia package:
```Julia

julia> ]
pkg> add MomentArithmetic

```
or for the most up to date version: 

```Julia
julia> ]
pkg> add https://github.com/AnderGray/MomentArithmetic.jl#master

```


## Usage

### Consistency Checking
```Julia
# Theoretical moment bounds from range
julia> A = Moments(mean = missing, var = missing, range = interval(2.3, 7))
moment: 	  ~ ( mean = [2.3,7.0], var = [0.0,5.5225] , range = [2.3,7.0] )

# Tighter variances when mean is given
julia> B = Moments(mean = 3, var = missing, range = interval(2.3, 7))
moment: 	  ~ ( mean  = 3, var = [0.0,2.8000000000000007] , range = [2.3,7.0] )

# Provided mean tightened with theoretical bound
julia> C = Moments(mean = interval(5, 10), var = missing, range = interval(2.3, 7))
moment: 	  ~ ( mean = [5.0,7.0], var = [0.0,5.4] , range = [2.3,7.0] )

# Gives error when provided information is outside theoretical bounds
julia> D = Moments(mean = interval(10, 30), var = missing, range = interval(2.3, 7))
ERROR: ArgumentError: Provided information not valid. Mean ∩ Range = ∅.
       [10, 30] ∩ [2.29999, 7] = ∅

# Variance outside theoretical bounds
julia> E = Moments(mean = 3, var = interval(15,18), range = interval(2.3, 7))
ERROR: ArgumentError: Provided information not valid. Variance ∩ VarBounds = ∅.
       [15, 18] ∩ [0, 2.80001] = ∅
```
### Arithmetic unary
```Julia

# Translation
julia> A + 2
moment: 	  ~ ( mean = [4.3,9.0], var = [0.0,5.5225] , range = [4.3,9.0] )

# Scaling
julia> A * 2
moment: 	  ~ ( mean = [4.6,14.0], var = [0.0,22.09] , range = [4.6,14.0] )

julia> 1/A
moment: 	  ~ ( mean = [0.14285714285714285,0.4347826086956522], var = [0.0,0.021305119401257674] , range = [0.14285714285714285,0.4347826086956522] )

julia> exp(A)
moment: 	  ~ ( mean = [9.974182454814718,1096.6331584284587], var = [0.0,295206.9325160221] , range = [9.974182454814718,1096.6331584284587] 

# Natural, also ln
julia> log(A)
moment: 	  ~ ( mean = [0.8329091229351039,1.9459101490553135], var = [0.0,0.3096928210361599] , range = [0.8329091229351039,1.9459101490553135] )

julia> log(A,10)
moment: 	  ~ ( mean = [0.3617278360175928,0.845098040014257], var = [0.0,0.058411688527944206] , range = [0.3617278360175928,0.845098040014257] )

julia> A^2
moment: 	  ~ ( mean = [5.289999999999999,49.0], var = [0.0,477.641025] , range = [5.289999999999999,49.0] )

julia> sqrt(A)
moment: 	  ~ ( mean = [1.51657508881031,2.6457513110645907], var = [0.0,0.3187597352261121] , range = [1.51657508881031,2.6457513110645907] )

julia> abs(A)
moment: 	  ~ ( mean = [2.3,7.0], var = [0.0,5.5225] , range = [2.3,7.0] )
```


### Arithmetic bivariate
```Julia
# Default is unknown interaction (Fréchet)
julia> A + B
moment: 	  ~ ( mean = [5.3,10.0], var = [0.0,16.187104249420315] , range = [4.6,14.0] )

julia> A - B
moment: 	  ~ ( mean = [-0.7000000000000002,4.0], var = [0.0,16.187104249420315] , range = [-4.7,4.7] )

julia> A * B
moment: 	  ~ ( mean = [5.289999999999999,24.93230212471016], var = [0.0,472.7449931126878] , range = [5.289999999999999,49.0] )

julia> A / B
moment: 	  ~ ( mean = [0.3285714285714285,3.043478260869566], var = [0.0,1.8426797770147763] , range = [0.3285714285714285,3.043478260869566] )

julia> min(A,B)
moment: 	  ~ ( mean = [2.3,3.0], var = [0.0,2.8000000000000007] , range = [2.3,7.0] )

julia> max(A,B)
moment: 	  ~ ( mean = [3.0,7.0], var = [0.0,5.5225] , range = [2.3,7.0] )

# Independence also possible
julia> sumIndep(A,B)
moment: 	  ~ ( mean = [5.3,10.0], var = [0.0,8.322500000000002] , range = [4.6,14.0] )

julia> subIndep(A,B)
moment: 	  ~ ( mean = [-0.7000000000000002,4.0], var = [0.0,8.322500000000002] , range = [-4.7,4.7] )

julia> multIndep(A,B)
moment: 	  ~ ( mean = [6.8999999999999995,21.0], var = [0.0,202.36550000000008] , range = [5.289999999999999,49.0] )

julia> divIndep(A,B)
moment: 	  ~ ( mean = [0.5590277777777776,3.009661835748793], var = [0.003501216042008588,1.8426797770147763] , range = [0.3285714285714285,3.043478260869566] )

julia> minIndep(A,B)
moment: 	  ~ ( mean = [2.3,3.0], var = [0.0,2.8000000000000007] , range = [2.3,7.0] )

julia> maxIndep(A,B)
moment: 	  ~ ( mean = [3.0,7.0], var = [0.0,5.5225] , range = [2.3,7.0] )

```

### Bounding p-box 

```Julia
julia> A = Moments(3, 1, interval(1,5))
moment: 	  ~ ( mean  = 3, var = 1 , range = [1.0,5.0] )

julia> p = makepbox(A);
julia> plot(p)       # may also do plot(A)

```

<img src="https://imgur.com/eqsdj7w.png" data-canonical-src="https://imgur.com/eqsdj7w.png" width="1500" />

The above p-box bounds all distribution functions with mean 3, variance 1 and range [1, 5]

## Funding

The authors would like to thank the gracious support from the EPSRC iCase studentship award 15220067. We also gratefully acknowledge funding from UKRI via the EPSRC and ESRC Centre for Doctoral Training in Risk and Uncertainty Quantification and Management in Complex Systems. This research was supported by the EPSRC through grant EP/R006768/1, which is acknowledgedfor its funding and support. This work has been carried out within the framework of the EUROfusion Consortium and has received funding from the Euratom research and training programme 2014-2018 and 2019-2020 under grant agreement No 633053. The views and opinions expressed herein do not necessarily reflect those of the European Commission 

(PDF) Distribution-free uncertainty propagation. Available from: https://www.researchgate.net/publication/352225779_Distribution-free_uncertainty_propagation
