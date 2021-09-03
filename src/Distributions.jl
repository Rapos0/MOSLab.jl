
abstract type DistributionFunc 
end

struct BoltzmanDist <: DistributionFunc
end

struct FermiDist <: DistributionFunc
end

struct BlakemoreDist <: DistributionFunc
    b::Number
end
BlakemoreDist() = BlakemoreDist(0.27)


function (b::BoltzmanDist)(E,Ef,T)
    return b(E-Ef,T)
end

function (b::BoltzmanDist)( ΔE,T)
    return b((ΔE)/(kb*T))
end

function (b::BoltzmanDist)( ΔE)
    return exp(ΔE)
end

function (b::FermiDist)(E,Ef,T)
    return b(E-Ef,T)
end

function (b::BlakemoreDist)(E,Ef,T)
    return b(E-Ef,T)
end

"""
# Arguments 
    eta 
Calculation One Half Order Fermi Integral of eta, using trapezoidal integral to approximat it
the accuracy was compared to Rational Chebychev approximations, and Gaussian Quadrature Calculations
with maximum absolute error smaller than 1e-12.
The defition of the Fermi Integral is such it approaches exp(-η) without normalization constants
"""
function trapFermi(eta)
    aj = 0.5
    range = 8.0
    if eta > 0.0
        range = sqrt(eta + 64.)
    end
    h = 0.5
    nmax = range / h
    sum = 0.
    for i = 1:nmax
        u = i * h
        ff = 2 * (u^(2 * aj + 1)) / (1 + exp(u * u - eta))
        sum = sum + ff
    end
    pol = 0.
    npole = 0
    bk1 = 0
    while bk1 <= 14 * pi
        npole = npole + 1
        bk = (2 * npole - 1) * pi
        rho = sqrt(eta * eta + bk * bk)
        t1 = 1
        t2 = 0
        if eta < 0
            tk = - aj * (atan(-bk / eta) + pi)
        elseif eta == 0
        tk = 0.5 * pi * aj
        else
            tk = aj * atan(bk / eta)
        end
        rk = - (rho^aj)
        tk = tk + 0.5 * atan(t2 / t1)
        if eta < 0
            rk = -rk
        end
        ak = (2. * pi / h) * sqrt(0.5 * (rho + eta))
        bk1 = (2. * pi / h) * sqrt(0.5 * (rho - eta))
        if bk1 <= (14. * pi)
            gama = exp(bk1)
            t1 = gama * sin(ak + tk) - sin(tk)
            t2 = 1.0 - 2.0 * gama * cos(ak) + gama * gama
            pol = pol + 4.0 * pi * rk * t1 / t2
        end 
    end 
    npole = npole - 1
    fdp = (sum * h + pol)*2/sqrt(pi)
end

"""
# Arguments 
    eta 
Calculation One Minus Half Order Fermi Integral of eta, using trapezoidal integral to approximat it
the accuracy was compared to Rational Chebychev approximations, and Gaussian Quadrature Calculations
with maximum absolute error smaller than 1e-12.
The defition of the Fermi Integral is such it approaches exp(-η) without normalization constants
"""
function DtrapFermi(eta)
    aj = -0.5
    range = 8.0
    if eta > 0.0
        range = sqrt(eta + 64.)
    end
    h = 0.5
    nmax = range / h
    sum = 0.
    for i = 1:nmax
        u = i * h
        ff = 2 * (u^(2 * aj + 1)) / (1 + exp(u * u - eta))
        sum = sum + ff
    end
    pol = 0.
    npole = 0
    bk1 = 0
    while bk1 <= 14 * pi
        npole = npole + 1
        bk = (2 * npole - 1) * pi
        rho = sqrt(eta * eta + bk * bk)
        t1 = 1
        t2 = 0
        if eta < 0
            tk = - aj * (atan(-bk / eta) + pi)
        elseif eta == 0
        tk = 0.5 * pi * aj
        else
            tk = aj * atan(bk / eta)
        end
        rk = - (rho^aj)
        tk = tk + 0.5 * atan(t2 / t1)
        if eta < 0
            rk = -rk
        end
        ak = (2. * pi / h) * sqrt(0.5 * (rho + eta))
        bk1 = (2. * pi / h) * sqrt(0.5 * (rho - eta))
        if bk1 <= (14. * pi)
            gama = exp(bk1)
            t1 = gama * sin(ak + tk) - sin(tk)
            t2 = 1.0 - 2.0 * gama * cos(ak) + gama * gama
            pol = pol + 4.0 * pi * rk * t1 / t2
        end 
    end 
    npole = npole - 1
    fdp = (sum * h + pol)*2/sqrt(pi)
end


"""
# Arguments 
    eta 
Calculation Three Half Order Fermi Integral of eta, using trapezoidal integral to approximat it
the accuracy was compared to Rational Chebychev approximations, and Gaussian Quadrature Calculations
with maximum absolute error smaller than 1e-12.
The defition of the Fermi Integral is such it approaches exp(-η) without normalization constants
"""
function trapFermi32(eta)
    aj = 1.5
    range = 8.0
    if eta > 0.0
        range = sqrt(eta + 64.)
    end
    h = 0.5
    nmax = range / h
    sum = 0.
    for i = 1:nmax
        u = i * h
        ff = 2 * (u^(2 * aj + 1)) / (1 + exp(u * u - eta))
        sum = sum + ff
    end
    pol = 0.
    npole = 0
    bk1 = 0
    while bk1 <= 14 * pi
        npole = npole + 1
        bk = (2 * npole - 1) * pi
        rho = sqrt(eta * eta + bk * bk)
        t1 = 1
        t2 = 0
        if eta < 0
            tk = - aj * (atan(-bk / eta) + pi)
        elseif eta == 0
            tk = 0.5 * pi * aj
        else
            tk = aj * atan(bk / eta)
        end
        rk = - (rho^aj)
        tk = tk + 0.5 * atan(t2 / t1)
        if eta < 0
            rk = -rk
        end
        ak = (2. * pi / h) * sqrt(0.5 * (rho + eta))
        bk1 = (2. * pi / h) * sqrt(0.5 * (rho - eta))
        if bk1 <= (14. * pi)
            gama = exp(bk1)
            t1 = gama * sin(ak + tk) - sin(tk)
            t2 = 1.0 - 2.0 * gama * cos(ak) + gama * gama
            pol = pol + 4.0 * pi * rk * t1 / t2
        end 
    end 
    npole = npole - 1
    fdp = (sum * h + pol)*2/sqrt(pi)
end

function Fermi12(x)
    Integrand(t,x) =  sqrt(t)/(exp(t-x)+1)
    t = Fun(t-> Integrand(t,x),0..10)
    return sum(t)
end


function (b::FermiDist)(ΔE,T)
    return b(ΔE/kb/T)
end
function (b::FermiDist)(ΔE)
    return trapFermi(ΔE)
end

function (b::BlakemoreDist)(ΔE,T)
    return b(ΔE/kb/T)
end
function (b::BlakemoreDist)(ΔE)
    return 1/(exp(-ΔE)+b.b)
end

function Primitive(b::BlakemoreDist,x)
    return log(b.b*exp(x)+1)/(b.b)
end

function Primitive(b::BoltzmanDist,x)
    return b(x)
end

function Primitive(b::FermiDist,x)
    return trapFermi32(x)
end