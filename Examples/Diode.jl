struct Diode
    Is
    n
    T
end

id()
Id(Va,Vc,m::Diode) = m.Is*(exp((Va-Vc)/kb/m.T/m.n)-1)
G(Va,Vc,m::Diode)