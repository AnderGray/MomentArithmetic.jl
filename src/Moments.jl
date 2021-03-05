######
# This file is part of the MomentArithmetic.jl package, for performing arithmetic operations
# between moments of uncertain numbers, when the moments and the uncertain numbers are only
# partially known.
#
#   Definition of Moment variable class
#
#           University of Liverpool, 
#           Institute for Risk and Unertainty
#
#                                           Author: Ander Gray, Scott Ferson
#                                           Email:  ander.gray@liverpool.ac.uk
#
#       To Do:
#
#               -> Mean var : cheb
#               -> Mean min: Markov, mean max: Markov
#               -> Mean min & max: Cantelli
#               -> Mean min & max var: Ferson
#   
#           (look in downs)
#   
#               -> Doesn't exits: bound and a variance
#
#               -> From moments, we can't get the range... unless we know something about the shapes
#               -> However possible option to truncate pba style. Give option. (Consitent throughout)
#
#               -> Univariate from Rowe  (found)
#
#               -> Arithmetic when dependency is known (?) Perf, op, correlation. Interval correlation coeff.
#               -> Dependeny. Calulate the covarience between Z and X. when z = x + y. Langewisch paper has some
#
#               
#
#               -> Units
#
######

import Base.+
import Base.-
import Base.*
import Base./
import Base.^

import Base.issubset, Base.intersect

using IntervalArithmetic, ProbabilityBoundsAnalysis

abstract type AbstractMoment end

include("Arithmetic.jl")
include("Utils.jl")


#global correlationDict = Dict()
#global unitDict = Dict()

#interval(x :: AbstractInterval, y :: AbstractInterval) = interval(x.lo, y.hi);

mutable struct Moments <: AbstractMoment

    mean ::  Union{Real, Interval}                   # Can also make it a pbox as an input
    var  ::  Union{Real, Interval}                   
    range :: Union{Missing, Interval}

    #CorreltedVars :: Union{Missing, Array{UInt64,1}} # For defining which variables this one is correlted too, for dependency tracking
    #Unit ::  Symbol   # Will add this later

    function Moments( mean = missing,  var = missing, range = missing)

        mean, var, range = make_consistent(mean, var, range)

        return new(mean, var, range)

    end
end

Moments(;mean = missing, var =missing, range =missing) = Moments(mean, var , range)

intervalM(mean = missing, var =missing, range =missing) =  Moments(mean, var , range)

intervalM(;mean =missing, var =missing, range =missing) =  Moments(mean, var , range)

function make_consistent(mean = missing,  var = missing, range = missing)


    if ismissing(mean);     mean    = interval(-Inf,Inf);   end
    if ismissing(var);      var     = interval(0,Inf);      end
    if ismissing(range);    range   = interval(-Inf,Inf);   end

    # mean constraint
    if (mean ∩ range == ∅); throw(ArgumentError("Provided information not valid. Mean ∩ Range = ∅.\n       $mean ∩ $range = ∅")); end
    if !( mean ⊆ range); 
        ml = max(left(mean), left(range));
        mh = min(right(mean), right(range));
        mean = interval(ml,mh);
    end

    # Var constraint
    maxVar = interval(0,Inf)
    if isfinite(range.hi) && isfinite(range.lo)
        maxRange = range.hi; minRange = range.lo; 
        meanMax = right(mean); meanMin = left(mean);

        v1 = (maxRange-minRange)*(maxRange-meanMin)-(maxRange-meanMin)*(maxRange-meanMin);
        v2 = (maxRange-minRange)*(maxRange-meanMax)-(maxRange-meanMax)*(maxRange-meanMax);
        v3 = 0; mid = (range.hi + range.lo)/2;

        if (mid ∈ mean); v3 =  (maxRange-minRange)*(maxRange-mid)-(maxRange-mid)*(maxRange-mid); end
        
        #println("v1: $v1");println("v2: $v2");println("v3: $v3")

        vh = max(v1, v2, v3); vl = 0;    
        maxVar = interval(vl,vh);
    end 

    if (var ∩ maxVar == ∅); throw(ArgumentError("Provided information not valid. Variance ∩ VarBounds = ∅.\n       $var ∩ $maxVar = ∅")); end
    if !(var ⊆ maxVar)
        vl = max(left(var), left(maxVar));
        vh = min(right(var), right(maxVar));
        var = interval(vl, vh)
    end

    return mean, var, range

end

function Base.show(io::IO, z::Moments)

    if (typeof(z.mean) <: AbstractInterval); statement1 = "mean = [$(z.mean.lo),$(z.mean.hi)]"; else statement1= "mean  = $(z.mean)";end
    if (typeof(z.var) <: AbstractInterval); statement2 = "var = [$(z.var.lo),$(z.var.hi)]"; else statement2 = "var = $(z.var)";end
    if ( !ismissing(z.range) ); statement3 = ", range = [$(z.range.lo),$(z.range.hi)]"; else; statement3 = "";end

    print(io, "moment: \t  ~ ( $statement1, $statement2 $statement3 )");

end




