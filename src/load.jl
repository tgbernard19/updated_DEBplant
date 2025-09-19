using Setfield
using Unitful, Plots, UnitfulRecipes, Microclimate, DataStructures, Flatten, FieldMetadata, 
      OrdinaryDiffEq, Photosynthesis, DynamicEnergyBudgets, Dates,
      DataStructures, ColorSchemes, JLD2, DimensionalData, FileIO

using DynamicEnergyBudgets: Records, PlottableRecords, scaling_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, parconv, w_V, build_vars, allometry_pars, 
      tempcorr, scaling_correction, photosynthesis

using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, g, mg, cm, m, s, hr, d, mol, mmol, μmol, σ

import Plots:px, pct, GridLayout

const STATELABELS = ["VS", "CS", "NS", "VR", "CR", "NR"]
const STATELABEL_SYMBOLS = Tuple(Symbol.(STATELABELS))
const STATE_INDEX = Dict(sym => idx for (idx, sym) in enumerate(STATELABEL_SYMBOLS))

state_index(label::Symbol) = STATE_INDEX[label]
state_index(label::AbstractString) = state_index(Symbol(label))

state_symbol(idx::Integer) = STATELABEL_SYMBOLS[idx]

function state_component(u, label::Symbol)
    if hasmethod(getindex, Tuple{typeof(u), Symbol})
        return u[label]
    elseif u isa AbstractArray
        return u[state_index(label)]
    elseif hasmethod(getproperty, Tuple{typeof(u), Symbol})
        return getproperty(u, label)
    else
        throw(ArgumentError("Unsupported state container $(typeof(u))"))
    end
end
state_component(u, label::AbstractString) = state_component(u, Symbol(label))

function state_set!(u, label::Symbol, value)
    if hasmethod(setindex!, Tuple{typeof(u), typeof(value), Symbol})
        u[label] = value
    elseif u isa AbstractArray
        u[state_index(label)] = value
    else
        throw(ArgumentError("Unsupported state container $(typeof(u))"))
    end
    return u
end
state_set!(u, label::AbstractString, value) = state_set!(u, Symbol(label), value)

function transect_from_saved(projectdir)
    locationspath = joinpath(projectdir, "data/locations.jld2")
    @load locationspath t1 t2 t3
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * u"hr"
    environments, tspan
end

function transect_from_netcdf(microclimdir::String, years, shade)
    envgrid = load_grid(microclimdir, years, shade);

    t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))
    t2 = MicroclimPoint(envgrid, CartesianIndex(60, 35))
    t3 = MicroclimPoint(envgrid, CartesianIndex(55, 35))

    save(joinpath(dir, "data/locations.jld2"), Dict("t1" => t1, "t2" => t2, "t3" => t3))

    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * u"hr"
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    environments, tspan
end

init_state(model::AbstractOrganism) = init_state(has_reserves.(define_organs(model, 1hr)), model)
init_state(::NTuple{2,HasCN}, model) = begin
    u = zeros(length(STATELABELS)) * mol
    state_set!(u, :VS, 0.2mg / (25.0g/mol))
    state_set!(u, :CS, 5.0mg  / (25.0g/mol))
    state_set!(u, :NS, 0.2mg  / (25.0g/mol))
    state_set!(u, :VR, 0.04mg / (25.0g/mol))
    state_set!(u, :CR, 1.0mg  / (25.0g/mol))
    state_set!(u, :NR, 0.04mg  / (25.0g/mol))
    u
end
init_state(::NTuple{2,HasCNE}, model) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.01]mol

# The zero crossing of allometry is the seed size.
set_allometry(model, state) = begin
    if :β0 in fieldnames(typeof(model.params[1].allometry_pars))
        model = @set model.params[1].allometry_pars.β0 = state_component(state, :VS) * w_V(model.shared) * 0.999999
    end
    if :β0 in fieldnames(typeof(model.params[2].allometry_pars))
        model = @set model.params[2].allometry_pars.β0 = state_component(state, :VR) * w_V(model.shared) * 0.999999
    end
    model
end

