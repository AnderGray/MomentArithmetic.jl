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


import Base: exp, sqrt, ^, log, ln

function env(x, y)
    MIN = min(left(x), left(y));
    MAX = max(right(x), right(y))
    return interval(MIN, MAX);
end

left(x :: AbstractMoment) = x.range.lo
right(x :: AbstractMoment) = x.range.hi

function rowe(x :: AbstractMoment, t :: Function)

    EX = x.mean; VX = x.var;
    LX = left(x); GX = right(x); 

    if typeof(EX) <: Interval;  throw(ArgumentError("Rowe mean must have precise moment inputs")); end    
    if typeof(VX) <: Interval;  throw(ArgumentError("Rowe mean must have precise moment inputs")); end    

    P = 1 / (1 + (EX - LX)^2/ VX)
    Q = 1 / (1 + (GX - EX)^2/ VX)

    MIN = P * t(LX) + (1 - P) * t(EX + VX/(EX - LX))
    MAX = Q * t(GX) + (1 - Q) * t(EX + VX/(EX - GX))

    return env(MIN, MAX)

end

function rowevar(x :: AbstractMoment, t :: Function, invT :: Function)
    
    EX = x.mean; VX = x.var;
    LX = left(x); GX = right(x); 

    if typeof(EX) <: Interval;  throw(ArgumentError("Rowe var must have precise moment inputs")); end    
    if typeof(VX) <: Interval;  throw(ArgumentError("Rowe var must have precise moment inputs")); end    

    v = invT(rowe(x, t))

    Lv = left(v); Gv = right(v);

    MIN = ((t(Lv) - t(LX))/(Lv - LX))^2 * (VX + (Lv - EX)^2)
    MAX = ((t(Gv) - t(GX))/(Gv - GX))^2 * (VX + (Gv - EX)^2)
    
    return env(MIN, MAX)

end


### 
#   Univariate transformations
###
-(x :: AbstractMoment) = Moments(-x.mean, x.var, -x.range)

function reciprocal(x :: AbstractMoment)
    f(x) = 1/x
    invF = *

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)

end

reciprocal(x :: Real) = 1/x

function exp(x :: AbstractMoment)

    MEAN = rowe(x, exp)
    VAR = rowevar(x, exp, log)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)

end

#^(x :: AbstractMoment) = throw(ArgumentError("rowe exponentiate must be implemented")) 

function ^(x :: AbstractMoment, a :: Real)

    if !(a == 2); throw(ArgumentError("Only rowe squared works so far")); end

    f(x) = x^2
    invF(x) = sqrt(x)

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
end

function sqrt(x :: AbstractMoment)

    f(x) = sqrt(x)
    invF(x) = x^2

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
end

function ln(x :: AbstractMoment)
    
    f(x) = log(x)
    invF(x) = exp(x)

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)
end

function log(x :: AbstractMoment, base = ℯ :: Real)
    if base == ℯ; ln(x);end

    f(x) = log(x, base)
    invF(x) = base^x

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
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
/(x :: Real, y :: AbstractMoment) = x * reciprocal(y);


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
    return Moments(zMean, zVar, zRange);
end

function subIndep(x :: AbstractMoment, y ::AbstractMoment)

    zMean  = x.mean  - y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range - y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);
end

function multIndep(x :: AbstractMoment, y ::AbstractMoment)

    zMean  = x.mean  * y.mean;
    zVar   = (x.mean^2 * y.var) + (y.mean^2 * x.var) + (x.var * y.var);
    zRange = x.range * y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);

end

function divIndep(x :: AbstractMoment, y ::AbstractMoment)
    multIndep(x,1/y);
end

+(x :: AbstractMoment, y ::AbstractMoment) = sumIndep(x,y);
-(x :: AbstractMoment, y ::AbstractMoment) = subIndep(x,y);
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