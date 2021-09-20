
struct SurrogateTransistor{T} <: TransistorModel
    W
    L
    T
    BaseModel::T
    Wrange
    Lrange
    Trange
    Voltagerange
    Idf
    gmf
    gdsf
end

function SurrogateTransistor(Wrange,Lrange,Trange,VoltageRange,model::TransistorModel)
    N_samples = 1000
    sampling_algorithm = LatinHypercubeSample()
    surrogate = LinearSurrogate
    opt = SRBF
    W = model.W
    L = model.L
    T = model.T0
    ## Create Ids Surrotage
    function fId(x)
        setW!(x[1],model)
        setL!(x[2],model)
        setT!(x[3],model)
        Id(x[4],x[5],x[6],model)
    end
    function fgm(x) 
        setW!(x[1],model)
        setL!(x[2],model)
        setT!(x[3],model)
        gm(x[4],x[5],x[6],model)
    end
    function fgds(x)
        setW!(x[1],model)
        setL!(x[2],model)
        setT!(x[3],model)
        gds(x[4],x[5],x[6],model)
    end
    lb = [Wrange[1],Lrange[1],Trange[1],VoltageRange[1],VoltageRange[1],VoltageRange[1]]
    ub = [Wrange[2],Lrange[2],Trange[2],VoltageRange[2],VoltageRange[2],VoltageRange[2]]
    sx = sample(N_samples,lb,ub,sampling_algorithm)

    sId = fId.(sx)
    suId = surrogate(sx,sId,lb,ub)
    surrogate_optimize(fId,opt(),lb,ub,suId,sampling_algorithm;maxiters=N_samples)

    sgm = fgm.(sx)
    sugm = surrogate(sx,sgm,lb,ub)
    surrogate_optimize(fgm,opt(),lb,ub,sugm,sampling_algorithm;maxiters=N_samples)

    sgds = fId.(sx)
    sugds = surrogate(sx,sgds,lb,ub)
    surrogate_optimize(fgds,opt(),lb,ub,sugds,sampling_algorithm;maxiters=N_samples)
    SurrogateTransistor(W,L,model.T0,model,Wrange,Lrange,Trange,VoltageRange,suId,sugm,sugds)

end

Id(Vg,Vd,Vs,mod::SurrogateTransistor) = mod.Idf([mod.W,mod.L,mod.T,Vg,Vd,Vs])
