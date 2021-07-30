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
#                                           Authors: Ander Gray, Scott Ferson
#                                           Email:  ander.gray@liverpool.ac.uk
######


##
#   Rowe mean
##

function roweFormula(EX, VX, LX , GX, t :: Function)

    P = 1 / (1 + (EX - LX)^2/ VX)
    Q = 1 / (1 + (GX - EX)^2/ VX)

    MIN = P * t(LX) + (1 - P) * t(EX + VX/(EX - LX))
    MAX = Q * t(GX) + (1 - Q) * t(EX + VX/(EX - GX))

    return env(MIN, MAX)

end


function rowe(x :: AbstractMoment, t :: Function; Nsub = Nsub[1])

    if !hasIntMoments(x); return roweNoSub(x,t); end

    EXs = x.mean; VXs = x.var;
    LX = left(x); GX = right(x);

    if typeof(EXs) <: Interval;  EXs = split(EXs, Nsub); end    # Sub-intervalise
    if typeof(VXs) <: Interval;  VXs = split(VXs, Nsub); end

    outs = [roweFormula(EX, VX, LX, GX, t) for EX in EXs, VX in VXs]

    return hull(outs[:])        # Same as env

end

# No sub-intervalisation
function roweNoSub(x :: AbstractMoment, t :: Function)

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

function rowevar(x :: AbstractMoment, t :: Function, invT :: Function; Nsub = Nsub[1])

    if !hasIntMoments(x); return rowevarNoSub(x,t, invT); end

    EXs = x.mean; VXs = x.var;
    LX = left(x); GX = right(x);

    if typeof(EXs) <: Interval;  EXs = split(EXs, Nsub); end    # Sub-intervalise
    if typeof(VXs) <: Interval;  VXs = split(VXs, Nsub); end

    outs = [rowevarFormula(EX, VX, LX, GX, t, invT) for EX in EXs, VX in VXs]

    out = hull(outs[:])
    if isnan(out); out = interval(0, Inf); end
    return out

end


# No sub-intervalisation
function rowevarNoSub(x :: AbstractMoment, t :: Function, invT :: Function)

    EX = x.mean; VX = x.var;
    LX = left(x); GX = right(x);

    return rowevarFormula(EX, VX, LX, GX, t, invT)

end

###
#   Univariate transformations
###
-(x :: AbstractMoment) = Moments(-x.mean, x.var, -x.range)

function reciprocal(x :: AbstractMoment)
    f(x) = 1/x
    invF(x) = 1/x

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)

end

reciprocal(x :: Real) = 1/x

function exp(x :: AbstractMoment)

    MEAN = rowe(x, exp)
    VAR = rowevar(x, exp, log)
    RANGE = exp(x.range)

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

function log(x = missing :: AbstractMoment, base = ℯ :: Real)
    if base == ℯ; return ln(x); end

    f(x) = log(base, x)
    invF(x) = base^x

    MEAN = rowe(x, f)
    VAR = rowevar(x, f, invF)
    RANGE = f(x.range)

    return Moments(MEAN, VAR, RANGE)
end

function abs(x :: AbstractMoment)
    if x.range > 0;
        EY = x.mean;
    elseif x.range < 0;
        EY = -x.mean;
    else
        sqvar = sqrt(x.var)
        EY1 = abs(x.mean)
        EY2 = EY1 + sqvar * (pi - atan(EY1/sqvar));
        EY = env(EY1, EY2)
    end

    Yvar = max(0, x.mean^2 + x.var - EY);
    Yrange = abs(x.range)

    return Moments(EY, Yvar, Yrange)


end

###
#   Scalar arithmetic
###
+(x :: AbstractMoment, y :: Real) = Moments(x.mean + y, x.var, x.range + y)
+(x :: Real, y:: AbstractMoment) = y + x;
-(x :: AbstractMoment, y :: Real) = x + (-y);
-(x :: Real, y :: AbstractMoment) = x + (-y);

*(x :: AbstractMoment, y :: Real) = Moments(x.mean * y, x.var * y^2, x.range * y)
*(x :: Real, y:: AbstractMoment) = y * x;
/(x :: AbstractMoment, y :: Real) = x * (1/y);
/(x :: Real, y :: AbstractMoment) = x * reciprocal(y);


###
#   Moment Arithmetic
###

###
# Independent
###
function sumIndep(x :: AbstractMoment, y :: AbstractMoment)

    zMean  = x.mean  + y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range + y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);
end

function subIndep(x :: AbstractMoment, y :: AbstractMoment)

    zMean  = x.mean  - y.mean;
    zVar   = x.var   + y.var;
    zRange = x.range - y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);
end

function multIndep(x :: AbstractMoment, y :: AbstractMoment)

    zMean  = x.mean  * y.mean;
    zVar   = (x.mean^2 * y.var) + (y.mean^2 * x.var) + (x.var * y.var);     # No subint required
    zRange = x.range * y.range;

    # Add dependent variables to dict
    return Moments(zMean, zVar, zRange);

end



function divIndep(x :: AbstractMoment, y :: AbstractMoment)
    multIndep(x,1/y);
end

function minIndep(x :: AbstractMoment, y :: AbstractMoment)

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

function maxIndep(x :: AbstractMoment, y :: AbstractMoment)

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

+(x :: AbstractMoment, y :: AbstractMoment) = sumFrechet(x,y);
-(x :: AbstractMoment, y :: AbstractMoment) = subFrechet(x,y);
*(x :: AbstractMoment, y :: AbstractMoment) = multFrechet(x,y);
/(x :: AbstractMoment, y :: AbstractMoment) = divFrechet(x,y);
min(x :: AbstractMoment, y :: AbstractMoment) = minFrechet(x,y);
max(x :: AbstractMoment, y :: AbstractMoment) = maxFrechet(x,y);

##
# With Intervals (Moments with just ranges)
##
+(x :: AbstractMoment, y :: Interval) = sumFrechet(x, Moments(missing, missing, y));
-(x :: AbstractMoment, y :: Interval) = subFrechet(x, Moments(missing, missing, y));
*(x :: AbstractMoment, y :: Interval) = multFrechet(x, Moments(missing, missing, y));
/(x :: AbstractMoment, y :: Interval) = divFrechet(x,  Moments(missing, missing, y));
min(x :: AbstractMoment, y :: Interval) = minFrechet(x,Moments(missing, missing, y));
max(x :: AbstractMoment, y :: Interval) = maxFrechet(x,Moments(missing, missing, y));

+(x :: Interval, y :: AbstractMoment) = sumFrechet(Moments(missing, missing, x), y);
-(x :: Interval, y :: AbstractMoment) = subFrechet(Moments(missing, missing, x), y);
*(x :: Interval, y :: AbstractMoment) = multFrechet(Moments(missing, missing, x),y);
/(x :: Interval, y :: AbstractMoment) = divFrechet(Moments(missing, missing, x), y);
min(x :: Interval, y :: AbstractMoment) = minFrechet(Moments(missing, missing, x),y);
max(x :: Interval, y :: AbstractMoment) = maxFrechet(Moments(missing, missing, x),y);

###
#   Dependence unknown
###
function sumFrechet(x :: AbstractMoment, y :: AbstractMoment)

    zMean = x.mean + y.mean

    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2

    zVarLb = max(zVarLb, 0)

    zVar = env(zVarLb, zVarUb);

    zRange = x.range + y.range;

    return Moments(zMean, zVar, zRange);
end

function subFrechet(x :: AbstractMoment, y :: AbstractMoment)

    zMean = x.mean - y.mean

    zVarLb = (sqrt(x.var) - sqrt(y.var))^2
    zVarUb = (sqrt(x.var) + sqrt(y.var))^2

    zVarLb = max(zVarLb, 0)

    zVar = env(zVarLb, zVarUb);

    zRange = x.range - y.range;

    return Moments(zMean, zVar, zRange);
end

# Requires alot of subintervaling for Gvar to have an effect
function multFrechet(x :: AbstractMoment, y :: AbstractMoment; Nsub = Nsub[1])

    EX = x.mean; EY = y.mean;
    VX = x.var; VY = y.var;

    zMeanLb = EX * EY - sqrt(VX * VY);
    zMeanUb = EX * EY + sqrt(VX * VY);
    zMean = env(zMeanLb, zMeanUb);

    if typeof(EX) <: Interval && !isvacuous(EX);  EX = split(EX, Nsub); end    # Sub-intervalise
    if typeof(EY) <: Interval && !isvacuous(EY);  EY = split(EY, Nsub); end

    zVar = [Gvar(Moments(ex ,VX ,x.range), Moments(ey, VY, y.range)) for ex in EX, ey in EY]

    zVar = hull(zVar[:])

    zRange = x.range * y.range;

    return Moments(zMean, zVar, zRange);
end

divFrechet(x :: AbstractMoment, y :: AbstractMoment) = multFrechet(x, 1/y)  #Not best possible


function minFrechet(x :: AbstractMoment, y :: AbstractMoment)

    zRange = min(x.range, y.range)
    zVar = env(0, max(x.var, y.var))

    EX = x.mean; EY = y.mean;

    LX = left(x.range); GX = right(x.range)
    LY = left(y.range); GY = right(y.range)

    zMean = Bmin(EX, LX, GX, EY, LY, GY);
    return Moments(zMean, zVar, zRange)

end

function maxFrechet(x :: AbstractMoment, y :: AbstractMoment)

    zRange = max(x.range, y.range)
    zVar = env(0, max(x.var, y.var))

    EX = x.mean; EY = y.mean;

    LX = left(x.range); GX = right(x.range)
    LY = left(y.range); GY = right(y.range)

    zMean = Bmax(EX, LX, GX, EY, LY, GY);
    return Moments(zMean, zVar, zRange)

end

###
#   Perfect Arithmetic
###

function sumPerfect(x :: Moments, y :: Moments)

    meanZ = x.mean + y.mean
    varZ  = (sqrt(x.var) + sqrt(y.var))^2
    rangeZ = x.range + y.range
    return Moments(meanZ, varZ, rangeZ)
end

function subPerfect(x :: Moments, y :: Moments)

    meanZ = x.mean - y.mean
    varZ  = (sqrt(x.var) - sqrt(y.var))^2
    rangeZ = x.range - y.range
    return Moments(meanZ, varZ, rangeZ)
end

###
#   Opposite Arithmetic
###

function sumOpposite(x :: Moments, y :: Moments)

    meanZ = x.mean + y.mean
    varZ  = (sqrt(x.var) - sqrt(y.var))^2
    rangeZ = x.range + y.range
    return Moments(meanZ, varZ, rangeZ)
end

function subOpposite(x :: Moments, y :: Moments)

    meanZ = x.mean - y.mean
    varZ  = (sqrt(x.var) + sqrt(y.var))^2
    rangeZ = x.range - y.range
    return Moments(meanZ, varZ, rangeZ)
end

###
#   known Dependence, or interval Dependence
###

function sumCor(x :: Moments, y :: Moments, corr :: Interval, nSub = 10)

    xVar = x.var; yVar = y.var; R = corr;

    if typeof(x.var) <: Interval; xVar = mince(xVar, nSub); end
    if typeof(y.var) <: Interval; yVar = mince(yVar, nSub); end
    #if typeof(corr) <: Interval; R = mince(R, nSub); end

    varZ = [xs + ys + 2 * sqrt(xs) * sqrt(ys) * corrs for xs in xVar, ys in yVar, corrs in R]

    meanZ = x.mean + y.mean
    rangeZ = x.range + y.range

    varZ = hull(varZ[:]);

    return Moments(meanZ, varZ, rangeZ)
end

sumCor(x :: Moments, y :: Moments, corr :: Real, nSub = 10) = sumCor(x,y, interval(corr), nSub)

function subCor(x :: Moments, y :: Moments, corr :: Interval, nSub = 10)

    xVar = x.var; yVar = y.var; R = corr;

    if typeof(x.var) <: Interval; xVar = mince(xVar, nSub); end
    if typeof(y.var) <: Interval; yVar = mince(yVar, nSub); end
    #if typeof(corr) <: Interval; R = mince(R, nSub); end

    varZ = [xs + ys - 2 * sqrt(xs) * sqrt(ys) * corrs for xs in xVar, ys in yVar, corrs in R]

    meanZ = x.mean - y.mean
    rangeZ = x.range - y.range

    varZ = hull(varZ[:]);

    return Moments(meanZ, varZ, rangeZ)
end

subCor(x :: Moments, y :: Moments, corr :: Real, nSub = 10) = subCor(x,y, interval(corr), nSub)

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

#   Looks like a complidated repeat variables problem, but only EX and EY need sub-ints
function Gvar(x :: AbstractMoment, y :: AbstractMoment)

    EX = x.mean; VX = y.var
    EY = y.mean; VY = y.var

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

    #VXY2 = EX^2*VX + EY^2*VY + 2*EX*EX2Y - 2*EY*EX2Y + 4*EX*EY^3 - 6*EX^2*EY^2 + 4*EX^3*EY - EXY^2 + 8*EX*EY*EXY - 4*EX^2*EXY - 4*EY^2*EXY -2*EX*EY*EX2 + EY^2*EX2 - 2*EX*EY*EY2 + EX^2*EY2 + 2*EY*EXY2 - 2*EX*EXY2 + EX2Y2

    return max(VXY, 0)

end

#=

    E11 = EXY - EX*EY
    E12 = EX2Y - EX2*EY + 2*EX^2*EY - 2*EX*EXY
    E21 = EXY2 - EX*EY2 + 2*EX*EY^2 - 2*EY*EXY
    E22 = -3*EX^2*EY^2 + EX2*EY^2 + EX^2*EY2 + 4*EX*EY*EXY - 2*EY*EX2Y - 2*EX*EXY2 + EX2Y2

    VXY = EX^2*VX + EY^2*VY + 2*EX*EY*E11 + 2*EX*E12 + 2*EY*E21 + E22 - E11^2

    A = EX; B = EY;
    C = EXY; D = EX2;
    E = EY2; F = EX2Y;
    G = EXY2; H = EX2Y2

    2 * A * B * (C - A*B) + 2 * A * (F -D*B +2 * A^2 * B - 2 * A * C) + 2 * B * (G - A * E + 2 * A * B^2 - 2 * B * C) + (-3 * A^2 * B^2 + D * B^2 + A^2 * E + 4*A*B*C -2*B*F - 2*A*G + H) - (C - A*B)^2

    simplify gives:

    2AF-2BF+4AB^3-6A^2B^2+4A^3B-C^2+8ABC-4A^2C-4B^2C-2ABD+B^2D-2ABE+A^2E+2BG-2AG+H

    2*EX*EX2Y - 2*EY*EX2Y + 4*EX*EY^3 - 6*EX^2*EY^2 + 4*EX^3*EY - EXY^2 + 8*EX*EY*EXY - 4*EX^2*EXY - 4*EY^2*EXY -2*EX*EY*EX2 + EY^2*EX2 - 2*EX*EY*EY2 + EX^2*EY2 + 2*EY*EXY2 - 2*EX*EXY2 + EX2Y2
=#

function multFrechetMean(x :: AbstractMoment, y :: AbstractMoment)
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
