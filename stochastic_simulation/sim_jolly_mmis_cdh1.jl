using ArgParse
using DifferentialEquations
using JLD2

parser = ArgParseSettings()
@add_arg_table parser begin
    "icfile"
        help = "input initial conditions"
    "psetfile"
        help = "input parameter set dictionary"
    "outfile"
        help = "output JLD2"
    "--cells"
        arg_type = Int
        default = 500
        help = "number of simulations"
    "--time"
        arg_type = Float64
        default = 500.
        help = "simulation endpoint"
    "--timechunk"
        arg_type = Float64
        default = 5.
        help = "simulation time per solve call"
end
args = parse_args(parser)

include("emt_jolly_mmis_cdh1.jl")

python_species_order = ["R_101", "X_SNAI", "X_SLUG", "R_200", "X_ZEB", "R_ZEB", "R_SNAI", "R_SLUG", "R_CDH", "C_ZEB_2", "C_ZEB_22", "C_ZEB_1", "C_ZEB_12", "C_ZEB_122", "C1_SLUG", "X_CDH"]
python_ics = include(args["icfile"])
ics = map(python_ics) do pyic
    ic = zeros(length(species_ids))
    for (s, conc) in enumerate(pyic)
        ic[species_ids[python_species_order[s]]] = conc
    end
    return ic
end

function filter_sols(sols)
    nvars, npts, nsols = size(sols)
    good_sols = []
    for i in 1:nsols
        if isassigned(sols[i][1,:], npts) == true
            push!(good_sols, i)
        end
    end
    return (good_sols, sols[:, :, good_sols])
end

function sim_pop(model, stoch, p, popsize, tmax, tstep=5.0, saveat=0.1)
    u0_blank = zeros(length(ics[1]))
    prob_sde = SDEProblem(model, stoch, u0_blank, (0., tstep), p)
    prob_func_cells(prob, i, repeat) = remake(prob, u0=ics[mod1(i, length(ics))])
    prob_cells = EnsembleProblem(prob_sde, prob_func=prob_func_cells)
    result = solve(prob_cells, saveat=saveat, maxiter=1e14, dt=0.0001, trajectories=popsize);
    _, result = filter_sols(result)
    last_t = tstep
    while last_t < tmax
        println("$(size(result)[3]) stable solutions at t = $last_t")
        succeeded = size(result)[3]
        prob_func_resume(prob, i, repeat) = remake(prob, u0=result[:, end, i])
        prob_resume = EnsembleProblem(prob_sde, prob_func=prob_func_resume)
        next_sol = solve(prob_resume, saveat=saveat, maxiter=1e14, dt=0.0001, trajectories=succeeded);
        good_sols, next_result = filter_sols(next_sol)
        result = hcat(result[:, :, good_sols], next_result[:, 2:end, :])
        last_t += tstep
    end
    println("$(size(result)[3]) stable solutions at final t = $last_t")
    result
end

noise = 0.2
p_array = copy(p_default)
p_array[90:(end - 1)] .= noise
pset = include(args["psetfile"])
for (p, value) in pset
    if p in keys(param_ids)
        p_array[param_ids[p]] = value
    end
end

println("Beginning simulation")
results = sim_pop(model, stoch, p_array, args["cells"], args["time"], args["timechunk"], 0.01)
save(args["outfile"], "results", results[1:16, 1:10:end, :])
