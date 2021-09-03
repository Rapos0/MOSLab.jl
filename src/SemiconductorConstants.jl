Dunit = 1.0
Eunit = 1.0
Tunit = 1.0
Funit = 1.0
Cunit = 1.0

kb = 8.6173303e-5*Eunit/Tunit
ϵ₀ = 8.85418782e-14*Funit/Dunit
q = 1.60217662e-19*Cunit

function defineConstants()
    global kb = 8.6173303e-5*Eunit/Tunit
    global ϵ₀ = 8.85418782e-14*Funit/Dunit
    global q = 1.60217662e-19*Cunit
end


function SemiUseUnits()
    global Dunit = u"cm"
    global Eunit = u"eV"
    global Tunit = u"K"
    global Funit = u"F"
    global Cunit = u"C"
    defineConstants()
end

function SemiUstrip()
    global Dunit = 1.0
    global Eunit = 1.0
    global Tunit = 1.0
    global Funit = 1.0
    global Cunit = 1.0
    defineConstants()
end
