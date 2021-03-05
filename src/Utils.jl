######
# This file is part of the MomentArithmetic.jl package, for performing arithmetic operations
# between moments of uncertain numbers, when the moments and the dependencies are only
# partially known.
#
#   Utility Functions
#
#           University of Liverpool, 
#           Institute for Risk and Unertainty
#
#                                           Author: Ander Gray, Scott Ferson
#                                           Email:  ander.gray@liverpool.ac.uk
######




issubset(x :: AbstractVector, y :: IntervalBox{N,T}) where {N,T} = ∈(x,y)
issubset(x :: IntervalBox, y :: Interval) = issubset(x[1],y)
issubset(x :: Interval, y :: IntervalBox) = issubset(x,y[1])

function intersect(x:: Real, y :: AbstractInterval)
    if x ∈ y
        return x
    end
    return ∅
end
intersect(x :: AbstractInterval, y ::Real) = intersect(y,x)


left(x :: Interval) = x.lo;         left(x :: Real) = x;
right(x :: Interval) = x.hi;        right(x :: Real) = x

hasIntMean(x :: Moments) = typeof(x.mean) <: Interval
hasIntVar(x :: Moments) = typeof(x.var) <: Interval
hasIntMoments(x :: Moments) = hasIntMean(x) || hasIntVar(x)
