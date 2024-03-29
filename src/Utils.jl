######
# This file is part of the MomentArithmetic.jl package, for performing arithmetic operations
# between moments of uncertain numbers, when the moments and the dependencies are only
# partially known.
#
#   Some usefull functions
#
#           University of Liverpool,
#           Institute for Risk and Unertainty
#
#                                           Authors: Ander Gray, Scott Ferson
#                                           Email:  ander.gray@liverpool.ac.uk
######


function env(x, y)
    MIN = min(left(x), left(y));
    MAX = max(right(x), right(y))
    return interval(MIN, MAX);
end


function env(x :: AbstractMoment, y :: AbstractMoment)
    return Moments(env(x.mean, y.mean), env(x.var, y.var), env(x.range, y.range));
end

left(x :: AbstractMoment) = x.range.lo
right(x :: AbstractMoment) = x.range.hi


function split(X :: Interval, splits :: Integer = n)
    ###
    #   Splits an interval into n subintervals. Julia's range function is very accurate.
    ###

    xRange = range(X.lo, stop=X.hi, length = splits+1)
    return interval.(xRange[1:end-1],xRange[2:end]);

end

#=
issubset(x :: AbstractVector, y :: IntervalBox{N,T}) where {N,T} = ∈(x,y)
issubset(x :: IntervalBox, y :: Interval) = issubset(x[1],y)
issubset(x :: Interval, y :: IntervalBox) = issubset(x,y[1])
=#
function intersect(x:: Real, y :: AbstractInterval)
    if x ∈ y
        return x
    end
    return ∅
end
intersect(x :: AbstractInterval, y ::Real) = intersect(y,x)


left(x :: Interval) = x.lo;         left(x :: Real) = x;
right(x :: Interval) = x.hi;        right(x :: Real) = x

hasIntMean(x :: AbstractMoment) = typeof(x.mean) <: Interval
hasIntVar(x :: AbstractMoment) = typeof(x.var) <: Interval
hasIntMoments(x :: AbstractMoment) = hasIntMean(x) || hasIntVar(x)

function dist(x :: Moments, y :: Moments)
    z = subPerfect(x, y)
    zAbs = abs(z)
    return zAbs.mean
end

area = dist

function areaWorst(x :: Moments, y :: Moments)
    z = subOpposite(x, y)
    zAbs = abs(z)
    return zAbs.mean
end

function areaFrechet(x :: Moments, y :: Moments)
    z = subFrechet(x, y)
    zAbs = abs(z)
    return zAbs.mean
end


function cov(x :: Moments, y :: Moments , cor :: Interval = interval(-1,1))
    if !(cor ⊆ interval(-1,1)); throw(ArgumentError("Correlation must be in [-1, 1]"));end
    return sqrt(x.var) * sqrt(y.var) * cor
end

cov(x, y, cor :: Real) = cov(x, y, interval(cor))
