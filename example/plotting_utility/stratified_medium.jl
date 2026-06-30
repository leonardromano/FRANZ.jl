############################################################################################################################################################
# function to run a model and return its output in a nice format for plotting

"Return the output of the model for a stratified medium"
function run_model_strat(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0, potential::Function=z->0.0)
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, cosθ=1.0, ϕ=0.0)

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 / 4π * input2func(model.E_51)

    # store output
    output = Dict{Symbol, Vector{Float64}}()
    output[:r]     = [pos[3] for pos in x[:, 1]]
    output[:v]     = [pos[3] for pos in v[:, 1]]
    output[:dA]    = [pos[3] for pos in n[:, 1]]
    output[:t]     = t[:, 1]
    output[:M]     = M[:, 1]
    output[:V]     = [x[i, 1]'n[i, 1]/3 for i in eachindex(x[:, 1])]
    output[:p]     = @. output[:M] * output[:v] / Num_SN(t[:, 1])
    output[:f_kin] = @. 0.5 * output[:M] * output[:v]^2 * km_s^2 * Msol / E_SN.(t[:, 1])
    output[:f_pot] = output[:M] * Msol .* potential.(output[:r]) ./ E_SN.(t[:, 1])

    # store the final outcome of the model
    # -1 -> NaN value encountered
    #  0 -> default (computation reached final time)
    #  1 -> singular surface normal
    #  2 -> z=0 reached
    output[:case] = f[:, 1]

    if output[:case][1] == 2
        # add correct state for z=0
        not_nan = .!(isnan.(output[:r]))
        N_end = length(output[:r][not_nan])

        dt = output[:r][N_end] / output[:v][N_end]

        # integrate
        output[:r][N_end+1]  = 1e-10
        output[:v][N_end+1]  = output[:v][N_end]
        output[:t][N_end+1]  = output[:t][N_end] + dt
        output[:dA][N_end+1] = output[:dA][N_end]
        output[:M][N_end+1]  = output[:M][N_end]
        output[:V][N_end+1]  = 1e-10 * output[:dA][N_end] / 3
        output[:p][N_end+1]  = output[:p][N_end]
        output[:f_kin][N_end+1]  = output[:f_kin][N_end]
        output[:f_pot][N_end+1]  = output[:M][N_end+1] * Msol * potential.(1e-10) / E_SN(output[:t][N_end+1])
    end

    return output
end

############################################################################################################################################################
# function for plotting the output

"Make plot showing radius, speed, mass, momentum per SN and energy efficiencies as a function of time for different models in a stratifed medium."
function create_figure_strat(models::Dict{String, Dict{Symbol, Vector{Float64}}}, model_names::Vector{String}, 
                             names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, tlim::Vector{<:Real})
    # create figure & axes
    fig = figure(figsize=(18, 10), num=1, clear=true)
    axs = fig.subplots(ncols=3, nrows=2)
    subplots_adjust(wspace=0.25, hspace=0.0)

    # reduce margins
    margins(0,0)

    # plot results
    for model in model_names
        axs[1,1].plot(models[model][:t], models[model][:r], color=c_model[model], linewidth=3, zorder=1)
        axs[1,2].plot(models[model][:t], models[model][:v], color=c_model[model], linewidth=3, zorder=1)
        axs[1,3].plot(models[model][:t], models[model][:M], color=c_model[model], linewidth=3, zorder=1)
        axs[2,1].plot(models[model][:t], models[model][:p], color=c_model[model], linewidth=3, zorder=1)
        axs[2,2].plot(models[model][:t], models[model][:f_kin], color=c_model[model], linewidth=3, zorder=1)
        axs[2,2].plot(models[model][:t], models[model][:f_pot], color=c_model[model], linestyle="--", linewidth=3, zorder=1)
        axs[2,3].plot([NaN], [NaN], color=c_model[model], label=" " * names[model], linewidth=3, zorder=1)
    end

    # Axis labels
    axs[1,1].set_ylabel(L"\mathrm{Shock \  Radius \ [pc]}", fontsize=fontsize)
    axs[1,2].set_ylabel(L"\mathrm{Shock \  Speed \ [km/s]}", fontsize=fontsize)
    axs[1,3].set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    axs[1,3].set_ylabel(L"\mathrm{Mass \ [M_{\odot} / sr]}", fontsize=fontsize)
    axs[2,1].set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    axs[2,1].set_ylabel(L"\mathrm{Momentum \ per \ SN \ [M_{\odot} \ km/s / sr]}", fontsize=fontsize)
    axs[2,2].set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    axs[2,2].set_ylabel(L"\mathrm{Energy \ Efficiency}", fontsize=fontsize)

    # add figure handles
    axs[1,1].text(0.05, 0.90, "(a)", transform=axs[1,1].transAxes, fontsize=fontsize)
    axs[1,2].text(0.05, 0.05, "(b)", transform=axs[1,2].transAxes, fontsize=fontsize)
    axs[1,3].text(0.05, 0.90, "(c)", transform=axs[1,3].transAxes, fontsize=fontsize)
    axs[2,1].text(0.05, 0.90, "(d)", transform=axs[2,1].transAxes, fontsize=fontsize)
    axs[2,2].text(0.05, 0.05, "(e)", transform=axs[2,2].transAxes, fontsize=fontsize)

    # decorate axes
    for ax in axs
        # set axis scales
        ax.set_xscale("log");ax.set_yscale("log")
        ax.set_xlim(tlim...)

        # size up tick parameters
        ax.tick_params(axis="both", labelsize=fontsize)
    end

    # remove x-tick labels in top row
    for i in 1:2
        axs[1, i].set_xticklabels([])
    end

    # energy fractions panel
    axs[2,2].plot([NaN], [NaN], color="black", label=" Kinetic Energy", linestyle="-", linewidth=3, zorder=1)
    axs[2,2].plot([NaN], [NaN], color="black", label=" Potential Energy", linestyle="--", linewidth=3, zorder=1)
    axs[2,2].legend(loc = "upper left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # 6th panel
    axs[2, 3].get_xaxis().set_visible(false)
    axs[2, 3].get_yaxis().set_visible(false)
    for side in ["top", "left", "right", "bottom"]
        axs[2, 3].spines[side].set_visible(false)
    end
    axs[2, 3].set_zorder(-1)
    axs[2, 3].legend(loc = "center left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, axs
end