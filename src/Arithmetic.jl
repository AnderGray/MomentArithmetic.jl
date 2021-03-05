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


import Base: exp, sqrt, ^, log, ln, +, -, *, /

global Nsub = 20

function env(x, y)
    MIN = min(left(x), left(y));
    MAX = max(right(x), right(y))
    return interval(MIN, MAX);
end

left(x :: Moments) = x.range.lo
right(x :: Moments) = x.range.hi


function Split(X :: Interval, splits :: Integer = n)
    ###
    #   Splits an interval into n subintervals. Julia's range function is very accurate.
    ###

    xRange = range(X.lo, stop=X.hi, length = splits+1)
    return interval.(xRange[1:end-1],xRange[2:end]);

end

##
#   Rowe mean
##

function roweFormula(EX, VX, LX , GX, t :: Function)

    P = 1 / (1 + (EX - LX)^2/ VX)
    Q = 1 / (1 + (GX .- EX)^2/ VX)

    MIN = P * t(LX) + (1 - P) * t(EX + VX/(EX - LX))
    MAX = Q * t(GX) + (1 - Q) * t(EX + VX/(EX - GX))

    return env(MIN, MAX)

end


function rowe(x :: Moments, t :: Function; Nsub = Nsub)

    if !hasIntMoments(x); return roweNoSub(x,t); end

    EXs = x.mean; VXs = x.var;
    LX = left(x); GX = right(x); 

    if typeof(EXs) <: Interval;  EXs = Split(EXs, Nsub); end    # Sub-intervalise
    if typeof(VXs) <: Interval;  VXs = Split(VXs, Nsub); end    

    outs = [roweFormula(EX, VX, LX, GX, t) for EX in EXs, VX in VXs]

    return hull(outs[:])        # Same as env

end

# No sub-intervalisation
function roweNoSub(x :: Moments, t :: Function)

    EX = x.mean; VX = x.var;
    LX = left(x); GX = right(x);

    return roweFormula(EX, VX, LX, GX, t)

end

###
#   Rowe Variance
##

function rowevarFormula(EX, VX, LX , GX, t :: Function, invT :: Function)


    v = invT(roweFormula(EX, VX, LX , GX, t))

    Lv = left(v); Gv = right(v);

    MIN = ((t(Lv) - t(LX))/(Lv - LX))^2 * (VX + (Lv - EX)^2)
    MAX = ((t(Gv) - t(GX))/(Gv - GX))^2 * (VX + (Gv - EX)^2)
    
    return env(MIN, MAX)

end

function rowevar(x :: Moments, t :: Function, invT :: Function; Nsub = Nsub)

    if !hasIntMoments(x); return rowevarNoSub(x,t, invT); end

    EXs = x.mean; VXs = x.var;
    LX = left(x); GX = right(x); 

    if typeof(EXs) <: Interval;  EXs = Split(EXs, Nsub); end    # Sub-intervalise
    if typeof(VXs) <: Interval;  VXs = Split(VXs, Nsub); end   

    outs = [rowevarFormula(EX, VX, LX, GX, t, invT) for EX in EXs, VX in VXs]
    
    return hull(outs[:])

end


# No sub-intervalisation
function rowevarNoSub(x :: Moments, t :: Function, invT :: Function)
    
    EX = x.mean; VX = x.var;
    LX = left(x); GX = right(x); 
    
    return rowevarFormula(EX, VX, LX, GX, t, invT)

end

### 
#   Univariate transformations
###
-(x :: Moments) = Moments(-x.mean, x.var, -x.range)

function reciprocal(x :: Moments)
    f(x) = 1/x
    invF(x) = 1/x

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)

end

reciprocal(x :: Real) = 1/x

function exp(x :: Moments)

    MEAN = rowe(x, exp)
    VAR = rowevar(x, exp, log)
    RANGE = exp(x.range)

    return Moments(MEAN, VAR, RANGE)

end

#^(x :: Moments) = throw(ArgumentError("rowe exponentiate must be implemented")) 

function ^(x :: Moments, a :: Real)

    if !(a == 2); throw(ArgumentError("Only rowe squared works so far")); end

    f(x) = x^2
    invF(x) = sqrt(x)

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
end

function sqrt(x :: Moments)

    f(x) = sqrt(x)
    invF(x) = x^2

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
end

function ln(x :: Moments)
    
    f(x) = log(x)
    invF(x) = exp(x)

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)
end

function log(x = missing :: Moments, base = ℯ :: Real)
    if base == ℯ; return ln(x); end

    f(x) = log(base, x)
    invF(x) = base^x

    MEAN = rowe(x, f) 
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)
    
    return Moments(MEAN, VAR, RANGE)
end


###
#   Scalar arithmetic 
###
+(x :: Moments, y :: Real) = Moments(x.mean + y, x.var, x.range + y)
+(x :: Real, y:: Moments) = y + x;
-(x :: Moments, y :: Real) = x + (-y);
-(x :: Real, y :: Moments) = x + (-y);

*(x :: Moments, y :: Real) = Moments(x.mean * y, x.var * y^2, x.range * y)
*(x :: Real, y:: Moments) = y * x;
/(x :: Moments, y :: Real) = x * (1/y);
/(x :: Real, y :: Moments) = x * reciprocal(y);


###
#   Moment Arithmetic
###

###
# Independent
###
function sumIndep(x :: Moments, y :: Moments)

    zMean  = x.mean  + y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range + y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);
end

function subIndep(x :: Moments, y :: Moments)

    zMean  = x.mean  - y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range - y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);
end

function multIndep(x :: Moments, y :: Moments)

    zMean  = x.mean  * y.mean;
    zVar   = (x.mean^2 * y.var) + (y.mean^2 * x.var) + (x.var * y.var);     # No subint required
    zRange = x.range * y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);

end



function divIndep(x :: Moments, y :: Moments)
    multIndep(x,1/y);
end

function minIndep(x :: Moments, y :: Moments)
    
    if x.range.hi < y.range.lo
        zMean = x.mean;
        zVar = x.var;
        zRange = min(x.range, y.range)
        return Moments(zMean, zVar, zRange)
    end

    if y.range.hi < x.range.lo
        zMean = y.mean;
        zVar = y.var;
        zRange = min(x.range, y.range)
        return Moments(zMean, zVar, zRange)
    end

    return minFrechet(x, y)
end

function maxIndep(x :: Moments, y :: Moments)
    
    if x.range.hi < y.range.lo
        zMean = y.mean;
        zVar = y.var;
        zRange = max(x.range, y.range)
        return Moments(zMean, zVar, zRange)
    end

    if y.range.hi < x.range.lo
        zMean = x.mean;
        zVar = x.var;
        zRange = max(x.range, y.range)
        return Moments(zMean, zVar, zRange)
    end

    return maxFrechet(x, y)
end

+(x :: Moments, y :: Moments) = sumIndep(x,y);
-(x :: Moments, y :: Moments) = subIndep(x,y);
*(x :: Moments, y :: Moments) = multIndep(x,y);
/(x :: Moments, y :: Moments) = divIndep(x,y);

##
# With Intervals (Moments with just ranges)
##
+(x :: Moments, y :: Interval) = sumIndep(x, Moment(missing, missing, y));
-(x :: Moments, y :: Interval) = subIndep(x, Moment(missing, missing, y));
*(x :: Moments, y :: Interval) = multIndep(x, Moment(missing, missing, y));
/(x :: Moments, y :: Interval) = divIndep(x,  Moment(missing, missing, y));

+(x :: Interval, y :: Moments) = sumIndep(Moment(missing, missing, x), y);
-(x :: Interval, y :: Moments) = subIndep(Moment(missing, missing, x), y);
*(x :: Interval, y :: Moments) = multIndep(Moment(missing, missing, x),y);
/(x :: Interval, y :: Moments) = divIndep(Moment(missing, missing, x), y);


###
#   Dependence unknown
###
function sumFrechet(x :: Moments, y :: Moments)

    zMean = x.mean + y.mean

    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2

    zVarLb = max(zVarLb, 0)

    zVar = env(zVarLb, zVarUb);

    zRange = x.range + y.range;

    return Moments(zMean, zVar, zRange);
end

function subFrechet(x :: Moments, y :: Moments)

    zMean = x.mean - y.mean

    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2

    zVarLb = max(zVarLb, 0)

    zVar = env(zVarLb, zVarUb);
    
    zRange = x.range - y.range;

    return Moments(zMean, zVar, zRange);
end

function multFrechet(x :: Moments, y :: Moments)

    zMeanLb = x.mean * y.mean - sqrt(x.var * y.var);
    zMeanUb = x.mean * y.mean + sqrt(x.var * y.var);
    zMean = env(zMeanLb, zMeanUb);
    
    zVar = Gvar(x, y)
    
    zRange = x.range * y.range;

    return Moments(zMean, zVar, zRange);
end

divFrechet(x :: Moments, y :: Moments) = multFrechet(x, 1/y)  #Not best possible


function minFrechet(x :: Moments, y :: Moments)

    zRange = min(x.range, y.range)
    zVar = env(0, max(x.var, y.var))

    EX = x.mean; EY = y.mean;

    LX = left(x.range); GX = right(x.range)
    LY = left(y.range); GY = right(y.range)

    zMean = Bmin(EX, LX, GX, EY, LY, GY);
    return Moments(zMean, zVar, zRange)

end

function maxFrechet(x :: Moments, y :: Moments)

    zRange = max(x.range, y.range)
    zVar = env(0, max(x.var, y.var))

    EX = x.mean; EY = y.mean;

    LX = left(x.range); GX = right(x.range)
    LY = left(y.range); GY = right(y.range)

    zMean = Bmax(EX, LX, GX, EY, LY, GY);
    return Moments(zMean, zVar, zRange)

end


###
#   known Dependence, or interval Dependence   
###


###
#   Bertsimas mean for min and max. Do not require sub-intervalisation due to env and intersect
###

function Bmin(EX, LX, GX, EY, LY, GY)

    B1 = env(min(EX, EY), min(LX, LY))
    B2 = EX + EY - env(max(EX, EY), max(GX, GY))

    return intersect(B1,B2)
end

function Bmax(EX, LX, GX, EY, LY, GY)

    B1 = env(max(EX, EY), max(GX, GY))
    B2 = EX + EY - env(min(EX, EY), min(LX, LY))
    
    return intersect(B1,B2)
end



###
#   Goodman formula for variance under unknown dependence
###

#   Looks like a complidated repeat variables problem, but onlt EX and EY need sub-ints
function Gvar(x :: Moments, y ::Moments)

    EX = x.mean; VX = y.mean
    EY = y.mean; VY = y.mean
    
    EXY = multFrechetMean(x, y)
    EX2Y = multFrechetMean(x^2, y)
    EXY2 = multFrechetMean(x, y^2)
    EX2Y2 = multFrechetMean(x^2, y^2)

    EX2 = (x^2).mean
    EY2 = (y^2).mean

    E11 = EXY - EX*EY
    E12 = EX2Y - EX2*EY + 2*EX^2*EY - 2*EX*EXY
    E21 = EXY2 - EX*EY2 + 2*EX*EY^2 - 2*EY*EXY
    E22 = -3*EX^2*EY^2 + EX2*EY^2 + EX^2*EY2 + 4*EX*EY*EXY - 2*EY*EX2Y - 2*EX*EXY2 + EX2Y2

    VXY = EX^2*VX + EY^2*VY + 2*EX*EY*E11 + 2*EX*E12 + 2*EY*E21 + E22 - E11^2

    return max(VXY,0)

end

function multFrechetMean(x :: Moments, y :: Moments)
    zMeanLb = x.mean * y.mean - sqrt(x.var * y.var);
    zMeanUb = x.mean * y.mean + sqrt(x.var * y.var);
    return env(zMeanLb, zMeanUb);
end



###
#   Woodpile
###

#=
###
#   Bertsimas with subintervalistion (does not require)
###
function minFrechetSub(x :: Moments, y ::Moments)

    zRange = min(x.range, y.range)
    zVar = env(0, max(x.var, y.var))

    EX = x.mean; EY = y.mean;

    LX = left(x.range); GX = right(x.range)
    LY = left(y.range); GY = right(y.range)

    if !hasIntMean(x) && !hasIntMean(y);
        zMean = Bmin(EX, LX, GX, EY, LY, GY);
        return Moments(zMean, zVar, zRange)
    end

    if typeof(EX) <: Interval;  EX = Split(EX, Nsub); end    
    if typeof(EY) <: Interval;  EY = Split(EY, Nsub); end   

    zMeans = [Bmin(EX, LX, GX, EY, LY, GY) for EX in EX, EY in EY];

    zMean = hull(zMeans[:])

    return Moments(zMean, zVar, zRange)

end



F(xmean, ymean, xvar, yvar) = (xmean^2 * yvar) + (ymean^2 * xvar) + (xvar * yvar)

xmean = interval(0,100)
xvar = interval(0,3)
ymean = interval(0,100)
yvar = interval(0,3)

xs = Split(xvar, 10)
ys = Split(yvar, 10)

outs  = [F(xmean, ymean, xvar, yvar) for xvar in xs, yvar in ys]
hull(outs[:])


=#