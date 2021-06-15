module MomentArithmetic

using IntervalArithmetic, ProbabilityBoundsAnalysis

import ProbabilityBoundsAnalysis: plot, makepbox
import IntervalArithmetic: issubset, intersect

using ProbabilityBoundsAnalysis: issubset

import Base.+
import Base.-
import Base.*
import Base./
import Base.^

import Base.issubset, Base.intersect
import Base: exp, sqrt, ^, log, +, -, *, /, abs, min, max

global Nsub = [20]        # Number of sub-intervals

abstract type AbstractMoment end

export

    # From Moments.jl
    Moments, intervalM,
    make_consistent,
    makepbox, plot,

    # From Arithmetic.jl
    rowe, roweNoSub, rowevar, rowevarNoSub,
    +, -, *, /, reciprocal, exp, ^, sqrt, ln, log,
    abs,
    sumIndep, subIndep, multIndep, divIndep, minIndep, maxIndep,
    sumFrechet, subFrechet, multFrechet, divFrechet, minFrechet,
    maxFrechet,

    sumPerfect, subPerfect, sumOpposite, subOpposite,

    # From Utils.jl
    env, left, right, split, issubset, intersect, hasIntMean, hasIntVar, hasIntMoments,

    area, dist,

    # From IntervalArithmetic.jl
    interval, Interval,

    # From ProbabilityBoundsAnalysis.jl
    pbox


include("Moments.jl")
include("Utils.jl")
include("Arithmetic.jl")

end
