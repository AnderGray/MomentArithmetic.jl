######
# This file is part of the MomentArithmetic.jl package, for performing arithmetic operations
# between moments of uncertain numbers, when the moments and the uncertain numbers are only
# partially known.
#
#   Definition of Moment Arithmetic
#
#           University of Liverpool, 
#           Institute for Risk and Unertainty
#
#                                           Author: Ander Gray, Scott Ferson
#                                           Email:  ander.gray@liverpool.ac.uk
######


### 
#   Univariate transformations
###
-(x :: AbstractMoment) = Moments(-x.mean, x.var, -x.range)

function reciprocal(x :: AbstractMoment)
    throw(ArgumentError("rowe reciprocate must be implemented")) 
end

function exp(x :: AbstractMoment)
    throw(ArgumentError("rowe exponentiate must be implemented")) 
end

^(x :: AbstractMoment) = throw(ArgumentError("rowe exponentiate must be implemented")) 

function ln(x :: AbstractMoment)
    throw(ArgumentError("rowe natural log must be implemented")) 
end

function log(x :: AbstractMoment, base :: Real)
    if base == ℯ; ln(x);end
    throw(ArgumentError("rowe log must be implemented")) 
end


###
#   Scalar arithmetic 
###
+(x :: AbstractMoment, y :: Real) = Moments(x.mean + y, x.var, x.range + y)
+(x :: Real, y::AbstractMoment) = y + x;
-(x :: AbstractMoment, y :: Real) = x + (-y);
-(x :: Real, y :: AbstractMoment) = x + (-y);

*(x :: AbstractMoment, y :: Real) = Moments(x.mean * y, x.var * y^2, x.range * y)
*(x :: Real, y::AbstractMoment) = y * x;
/(x :: AbstractMoment, y :: Real) = x * (1/y);
/(x :: Real, y::AbstractMoment) = y * reciprocal(x);


###
#   Moment Arithmetic
###

###
# Independent
###
function sumIndep(x :: AbstractMoment, y ::AbstractMoment)

    zMean  = x.mean  + y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range + y.range;

    # Add dependent variables to dict
    return Moments(zMean,zVar,zRange);
end

function subIndep(x :: AbstractMoment, y ::AbstractMoment)

    zMean  = x.mean  - y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range - y.range;

    # Add dependent variables to dict
    return Moments(zMean,zVar,zRange);
end

function multIndep(x :: AbstractMoment, y ::AbstractMoment)

    zMean  = x.mean  * y.mean;
    zVar   = (x.mean^2 * y.var) + (y.mean^2 * x.var) + (x.var * y.var);
    zRange = x.range * y.range;

    # Add dependent variables to dict
    return Moments(zMean,zVar,zRange);

end

function divIndep(x :: AbstractMoment, y ::AbstractMoment)
    multIndep(x,1/y);
end

+(x :: AbstractMoment, y ::AbstractMoment) = sumIndep(x,y);
-(x :: AbstractMoment, y ::AbstractMoment) = subIndep(x.y);
*(x :: AbstractMoment, y ::AbstractMoment) = multIndep(x,y);
/(x :: AbstractMoment, y ::AbstractMoment) = divIndep(x,y);

###
#   Dependence unknown
###
function sumFrechet(x :: AbstractMoment, y ::AbstractMoment)

    zMean = x.mean + y.mean
    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2

    if typeof(zVarLb) <: AbstractInterval; zVarLb = zVarLb.lo; end
    if typeof(zVarUb) <: AbstractInterval; zVarUb = zVarUb.hi; end

    zVar = interval(zVarLb, zVarUb);
    
    zRange = x.range + y.range;

    return Moments(zMean, zVar, zRange);
end

function subFrechet(x :: AbstractMoment, y ::AbstractMoment)

    zMean = x.mean - y.mean

    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2
    zVar = interval(zVarLb, zVarUb);
    
    zRange = x.range - y.range;

    return Moments(zMean, zVar, zRange);
end

function multFrechet(x :: AbstractMoment, y ::AbstractMoment)

    zMeanLb = x.mean * y.mean - sqrt(x.var * y.var);
    zMeanUb = x.mean * y.mean + sqrt(x.var * y.var);
    zMean = interval(zMeanLb, zMeanUb);
    
    # For Var
    throw(ArgumentError("Goodman variance must be implemented")) 
    
    zVar = interval(zVarLb, zVarUb);
    
    zRange = x.range * y.range;

    return Moments(zMean, zVar, zRange);
end

divFrechet(x :: AbstractMoment, y ::AbstractMoment) = multFrechet(x, 1/y)  #Not best possible


###
#   known Dependence, or interval Dependence   
###