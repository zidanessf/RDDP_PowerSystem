using PowerModels,PrettyPrinting,JuMP
silence()
case = PowerModels.parse_file(joinpath(@__DIR__, "../input/case300.m"))
global n = 0
for gen in keys(case["gen"])
    if length(case["gen"][gen]["cost"]) == 0
        global n = n + 1
    end
end
println(n)
