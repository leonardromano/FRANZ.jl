############################################################################################################################################################
# functions to run a model and return its output in a nice format for plotting

"Return the output of the model for a stratified medium"
function run_model_single_fil(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0, potential::Function=x->0.0)
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, cosθ=0.0, ϕ=0.0)

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 / 4π * input2func(model.E_51)

    # store output
    output = Dict{Symbol, Vector{Float64}}()
    output[:t]     = t[:, 1]
    output[:r]     = [pos[1] for pos in x[:, 1]]
    output[:v]     = [pos[1] for pos in v[:, 1]]
    output[:dA]    = [pos[1] for pos in n[:, 1]]
    output[:M]     = M[:, 1]
    output[:V]     = [x[i, 1]'n[i, 1]/3 for i in eachindex(x[:, 1])]
    output[:n_H]   = output[:M] * Msol ./ (output[:V] * pc^3) / μ 
    output[:p]     = @. output[:M] * output[:v] / Num_SN(t[:, 1])
    output[:f_kin] = @. 0.5 * output[:M] * output[:v]^2 * km_s^2 * Msol / E_SN.(t[:, 1])
    output[:f_pot] = output[:M] * Msol .* potential.(x[:, 1]) ./ E_SN.(t[:, 1])

    return output
end

"Return the output of the model for a slice through a 3D medium"
function run_model_filaments_slice(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0,
    cosθ::Union{Real, AbstractRange{<:Real}, Vector{<:Real}}=LinRange(-1, 1, 5), ϕ::Union{Real, AbstractRange{<:Real}, Vector{<:Real}}=LinRange(0, 2π, 5))
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, cosθ=cosθ, ϕ=ϕ)

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 / 4π * input2func(model.E_51)

    output = Dict{Symbol, Any}()
    output[:t] = t
    output[:x] = x
    output[:n] = n
    output[:v] = v
    output[:p] = M .* v ./ Num_SN.(t)
    output[:case] = f

    return output
end

"Return the output of the model for a 3D structured medium."
function run_model_filaments(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0, potential::Function=x->0.0, Nside::Integer=4)
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, map=HealpixMap{Float64, NestedOrder}(Nside))

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 * input2func(model.E_51)

    # replace NaN values with more reasonable values
    for idir in eachindex(f)
        if f[idir] == 1
            # post collapse
            not_nan = .!(isnan.(t[:, idir]))
            iend = length(t[not_nan, idir])

            # assign reasonable values (frozen in place, zero volume)
            t[(iend+1):end, idir] = time[(iend+1):end]  
            @views fill!(x[(iend+1):end, idir], x[iend, idir])
            @views fill!(v[(iend+1):end, idir], zeros(3))
            @views fill!(n[(iend+1):end, idir], normalize(x[iend, idir]) * 1e-10)
            M[(iend+1):end, idir] .= M[iend, idir]
        elseif f[idir] == 2
            # post fountain-flow
            not_nan = .!(isnan.(t[:, idir]))
            iend = length(t[not_nan, idir])

            # assign reasonable values (frozen in place, non-zero volume)
            t[(iend+1):end, idir] = time[(iend+1):end]
            @views fill!(x[(iend+1):end, idir], x[iend, idir])
            @views fill!(v[(iend+1):end, idir], zeros(3))
            @views fill!(n[(iend+1):end, idir], n[iend, idir])
            M[(iend+1):end, idir] .= M[iend, idir]
        end
    end

    # Number of snapshots
    Ntime = length(time)
    dΩ = 4π / (Nside^2 * 12) # area associated with Healpix NSIDE

    # geometry
    r_shear     = zeros(Ntime, 3)
    r_vol       = zeros(Ntime, 3)
    minor_ratio = zeros(Ntime)
    major_ratio = zeros(Ntime)
    ϕ_minor     = zeros(Ntime)
    cosθ_minor  = zeros(Ntime)
    ϕ_major     = zeros(Ntime)
    cosθ_major  = zeros(Ntime)

    # dynamics
    vel_shear       = zeros(Ntime)
    v_t_shear       = zeros(Ntime)
    Mass_shear      = zeros(Ntime)
    Volume_shear    = zeros(Ntime)
    Density_shear   = zeros(Ntime)
    Area_shear      = zeros(Ntime)
    p_shear         = zeros(Ntime)
    p_t             = zeros(Ntime)
    f_kin_rad_shear = zeros(Ntime)
    f_kin_t_shear   = zeros(Ntime)
    f_pot_shear     = zeros(Ntime)

    # get center and geometry ellipsoid
    x_s, I_s = compute_deformation(x)

    # get geometry
    for it in eachindex(time)
        abc, dirs = get_ellipsoid(I_s[it, :, :])

        # explicitly normalize eigenvectors
        for dir in dirs 
            normalize!(dir)
        end

        # get cylindrical coordinates
        z = [xi[3] for xi in x[it, :]]
        R = [norm(xi[1:2]) for xi in x[it, :]]

        # get SNR radius
        r_shear[it, 1] = prod(abc)^(1/3)
        r_shear[it, 2] = minimum(abc)
        r_shear[it, 3] = maximum(abc)

        minor_ratio[it] = abc[3] / abc[1]
        major_ratio[it] = abc[2] / abc[1]

        # get angles
        e_φ = [0, 1, 0]
        e_R = [1, 0, 0]

        # get minor axis
        cosθ    = dirs[3][3]
        φ       = atan(e_R'dirs[3], e_φ'dirs[3]) * 180/π
        ϕ_minor[it]    = forward_frame(φ, acosd(cosθ))[1]
        cosθ_minor[it] = cosθ

        # get major axis
        cosθ    = dirs[3][1]
        φ       = atan(e_R'dirs[1], e_φ'dirs[1]) * 180/π
        ϕ_major[it]    = forward_frame(φ, acosd(cosθ))[1]
        cosθ_major[it] = cosθ

        # get Mass and volume
        Mass_shear[it]   = sum(M[it, :]) * dΩ
        Volume_shear[it] = sum([dot(x[it, ihp], n[it, ihp]) / 3 for ihp in eachindex(x[it, :])]) * dΩ
        Density_shear[it] = Mass_shear[it] * Msol / Volume_shear[it] / pc^3 / μ 
        Area_shear[it]   = sum([norm(dA) for dA in n[it, :]]) * dΩ
        r_vol[it, 1]     = (3 * Volume_shear[it] / 4π)^(1/3)
        r_c              = [norm(xi) for xi in x[it, :]]
        r_vol[it, 2]     = minimum(r_c) 
        r_vol[it, 3]     = maximum(r_c) 

        # get velocities
        e_rad = [normalize(dA) for dA in n[it, :]]

        # get radial and transversal velocities for all directions
        v_rad = [e_rad[ihp]'v[it, ihp] for ihp in eachindex(e_rad)]
        v_t   = [norm(v[it, ihp] - v_rad[ihp] * e_rad[ihp]) for ihp in eachindex(v_rad)]

        # compute radial and transversal momentum
        p_rad  = sum(M[it, :] .* v_rad) * dΩ
        p_perp = sum(M[it, :] .* v_t) * dΩ

        # get mass-weighted mean velocities
        vel_shear[it] = p_rad / Mass_shear[it]
        v_t_shear[it] = p_perp / Mass_shear[it]

        # get Momentum
        p_shear[it] = p_rad / Num_SN(time[it])
        p_t[it]     = p_perp / Num_SN(time[it])

        # get kinetic energy
        f_kin_rad_shear[it] = 0.5 * sum(M[it, :] .* v_rad.^2) * dΩ * km_s^2 * Msol / E_SN(time[it])
        f_kin_t_shear[it]   = 0.5 * sum(M[it, :] .* v_t.^2) * dΩ * km_s^2 * Msol / E_SN(time[it])

        # get potential energy
        f_pot_shear[it] = sum(M[it, :] .* potential.(x[it, :])) * Msol * dΩ / E_SN(time[it])
    end

    output = Dict{Symbol, Array{Float64}}()
    # geometry
    output[:r]           = r_shear
    output[:minor_ratio] = minor_ratio
    output[:major_ratio] = major_ratio
    output[:ϕ_minor]     = ϕ_minor     
    output[:cosθ_minor]  = cosθ_minor 
    output[:ϕ_major]     = ϕ_major     
    output[:cosθ_major]  = cosθ_major
    output[:r_vol]       = r_vol

    # Dynamics
    output[:v]       = vel_shear
    output[:v_t]     = v_t_shear
    output[:M]       = Mass_shear
    output[:V]       = Volume_shear
    output[:n_H]     = Density_shear
    output[:A]       = Area_shear
    output[:p]       = p_shear
    output[:p_t]     = p_t
    output[:f_kin]   = f_kin_rad_shear
    output[:f_kin_t] = f_kin_t_shear
    output[:f_pot]   = f_pot_shear

    return output
end

############################################################################################################################################################
# functions for setting up the environment

"Return the distance vector to the filament"
function get_distance_from_filament(x::Vector{<:Real}, x_fil::Vector{<:Real}, dir_fil::Vector{<:Real})
    dx = x_fil - x
    e_fil = normalize(dir_fil)
    
    return dx - dx'e_fil * e_fil
end

"Return the density profile of an isothermal infinite filament."
function ρ_filament(x::Vector{<:Real}; R_fil::Real, x_fil::Vector{<:Real}, dir_fil::Vector{<:Real})
    # make sure direction vector is normalized
    normalize!(dir_fil)

    # get the distance vector and its squared norm
    dx = get_distance_from_filament(x, x_fil, dir_fil) / R_fil
    r2 = dx'dx

    return 1 / (1 + r2)^2
end

"Return the gravitational acceleration from an isothermal infinite filament."
function g_filament(x::Vector{<:Real}; R_fil::Real, x_fil::Vector{<:Real}, dir_fil::Vector{<:Real})
    # make sure direction vector is normalized
    normalize!(dir_fil)

    # get the distance vector and its squared norm
    dx = get_distance_from_filament(x, x_fil, dir_fil) / R_fil
    r2 = dx'dx

    return 4 * dx / (1 + r2)
end

"Return the gravitational potential of an isothermal infinite filament."
function Φ_filament(x::Vector{<:Real}; R_fil::Real, x_fil::Vector{<:Real}, dir_fil::Vector{<:Real})
    # make sure direction vector is normalized
    normalize!(dir_fil)
    
    # get the distance vector and its squared norm
    dx = get_distance_from_filament(x, x_fil, dir_fil) / R_fil
    r2 = dx'dx

    return 2 * log(1 + r2)
end

############################################################################################################################################################
# functions for plotting the output

"Make plot showing radius, speed, mass, momentum per SN and energy efficiencies as a function of time for different models in a stratifed medium."
function create_figure_single_fil(models::Dict{String, Dict{Symbol, Vector{Float64}}}, model_names::Vector{String}, 
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
        axs[1,3].plot(models[model][:t], models[model][:M], color=c_model[model], linewidth=3, zorder=1, label=" " * names[model])
        axs[2,1].plot(models[model][:t], models[model][:p], color=c_model[model], linewidth=3, zorder=1)
        axs[2,2].plot(models[model][:t], models[model][:f_kin], color=c_model[model], linewidth=3, zorder=1)
        axs[2,2].plot(models[model][:t], models[model][:f_pot], color=c_model[model], linestyle="--", linewidth=3, zorder=1)
        axs[2,3].plot(models[model][:t], models[model][:n_H], color=c_model[model], linewidth=3, zorder=1)
    end

    # y-axis labels
    axs[1,1].set_ylabel(L"\mathrm{Shock \  Radius \ [pc]}", fontsize=fontsize)
    axs[1,2].set_ylabel(L"\mathrm{Shock \  Speed \ [km/s]}", fontsize=fontsize)
    axs[1,3].set_ylabel(L"\mathrm{Mass \ [M_{\odot} / sr]}", fontsize=fontsize)
    axs[2,1].set_ylabel(L"\mathrm{Momentum \ per \ SN \ [M_{\odot}/ sr \ km/s]}", fontsize=fontsize)  
    axs[2,2].set_ylabel(L"\mathrm{Energy \ Efficiency}", fontsize=fontsize)
    axs[2,3].set_ylabel(L"\mathrm{Density \ [cm^{-3}]}", fontsize=fontsize)

    # add figure handles
    axs[1,1].text(0.05, 0.05, "(a)", transform=axs[1,1].transAxes, fontsize=fontsize)
    axs[1,2].text(0.90, 0.90, "(b)", transform=axs[1,2].transAxes, fontsize=fontsize)
    axs[1,3].text(0.05, 0.90, "(c)", transform=axs[1,3].transAxes, fontsize=fontsize)
    axs[2,1].text(0.90, 0.05, "(d)", transform=axs[2,1].transAxes, fontsize=fontsize)
    axs[2,2].text(0.05, 0.05, "(e)", transform=axs[2,2].transAxes, fontsize=fontsize)
    axs[2,3].text(0.05, 0.90, "(f)", transform=axs[2,3].transAxes, fontsize=fontsize)

    # decorate axes
    for ax in axs
        # set axis scales
        ax.set_xscale("log");ax.set_yscale("log")
        ax.set_xlim(tlim...)

        # size up tick parameters
        ax.tick_params(axis="both", labelsize=fontsize)
    end

    # remove x-tick labels in top row
    for i in 1:3
        # no x-tick labels for the top row
        axs[1, i].set_xticklabels([])

        # axis label for bottom row
        axs[2,i].set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    end

    # energy fractions panel
    axs[2,2].plot([NaN], [NaN], color="black", label=" Kinetic", linestyle="-", linewidth=3, zorder=1)
    axs[2,2].plot([NaN], [NaN], color="black", label=" Potential", linestyle="--", linewidth=3, zorder=1)
    axs[2,2].legend(loc = "upper right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # legend items
    axs[1,3].plot([NaN], [NaN], color="black", label=" Uniform Medium", linestyle=":", linewidth=3, zorder=1)
    axs[1,3].legend(loc = "lower right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, axs
end

"Plot showing the contour lines of the shock surface at different times for different models. Also show velocity vectors as arrows."
function plot_contours_filaments(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Any}}, model_names::Vector{String}, 
                                 names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, times_frame::Vector{<:Real}, x_fil::Vector{<:Real}; 
                                 L_frame::Vector{<:Real}=[600, 600, 600], vel_scale::Vector{<:Real}=[250, 100, 75], N_arrow_sampling::Integer=4)

    # create figure
    fig = figure(figsize=(18, 12), num=1, clear=true)
    axs = fig.subplots(ncols=3, nrows=2)
    subplots_adjust(wspace=0.01, hspace=0.01)

    # for which subset to plot arrows?
    Nφ        = size(models["SN_xy"][:x])[2]
    Δi_sample = max(div(Nφ, N_arrow_sampling), 1)
    samples   = 1:Δi_sample:Nφ

    # get orbital times for each panel
    i_orb = get_timepoints(time, times_frame)

    for i_time in eachindex(times_frame)
        # get axis and center
        i_out = i_orb[i_time]

        # top row: "xy-plane"
        ax = axs[1, i_time]
        
        for model in model_names
            mdl = model * "_xy" 
            
            ax.plot([pos[1] for pos in models[mdl][:x][i_out, :]], [pos[2] for pos in models[mdl][:x][i_out, :]], color=c_model[model], linewidth=3)
            ax.quiver([pos[1] for pos in models[mdl][:x][i_out, samples]], [pos[2] for pos in models[mdl][:x][i_out, samples]], 
            [pos[1] for pos in models[mdl][:v][i_out, samples]], [pos[2] for pos in models[mdl][:v][i_out, samples]], width=0.005, color=c_model[model], scale=vel_scale[i_time])
        end

        ax.set_xlim(-0.5L_frame[i_time], 0.5L_frame[i_time])
        ax.set_ylim(-0.5L_frame[i_time], 0.5L_frame[i_time])

        # decorate axes
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.annotate(L"\mathrm{t = %$(times_frame[i_time]) \ Myr}", xy=[-0.45, 0.425] * L_frame[i_time], fontsize=fontsize, color="black")

        x_ruler = [-0.45, -0.45 + 1/6] * L_frame[i_time]; y_ruler = [-0.45, -0.45] * L_frame[i_time]
        ax.plot(x_ruler, y_ruler, color="black", linewidth=2)
        ax.annotate("$(round(Int, 1/6 * L_frame[i_time])) pc", xy=(x_ruler[1]-0.005 * L_frame[i_time], mean(y_ruler)+0.025 * L_frame[i_time]), fontsize=fontsize, color="black")  

        # bottom row: "xz-plane"
        ax = axs[2, i_time]
        
        for model in model_names
            mdl = model * "_xz" 
            
            ax.plot([pos[1] for pos in models[mdl][:x][i_out, :]], [pos[3] for pos in models[mdl][:x][i_out, :]], color=c_model[model], linewidth=3)
            ax.quiver([pos[1] for pos in models[mdl][:x][i_out, samples]], [pos[3] for pos in models[mdl][:x][i_out, samples]], 
            [pos[1] for pos in models[mdl][:v][i_out, samples]], [pos[3] for pos in models[mdl][:v][i_out, samples]], width=0.005, color=c_model[model], scale=vel_scale[i_time])
        end

        ax.set_xlim(-0.005L_frame[i_time], 0.5L_frame[i_time])
        ax.set_ylim(-0.005L_frame[i_time], 0.5L_frame[i_time])

        # decorate axes
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.annotate(L"\mathrm{t = %$(times_frame[i_time]) \ Myr}", xy=[0.025, 0.45] * L_frame[i_time], fontsize=fontsize, color="black")

        x_ruler = [0.375, 0.375 + 1/12] * L_frame[i_time]; y_ruler = [0.025, 0.025] * L_frame[i_time]
        ax.plot(x_ruler, y_ruler, color="black", linewidth=2)
        ax.annotate("$(round(Int, 1/12 * L_frame[i_time])) pc", xy=(x_ruler[1]+0.005 * L_frame[i_time], mean(y_ruler)+0.0125 * L_frame[i_time]), fontsize=fontsize, color="black")  
    end

    # add legend
    axs[2,1].plot([NaN], [NaN], color=c_model["SN"], label=" " * names["SN"], linestyle="-", linewidth=3, zorder=1)
    axs[2,1].plot([NaN], [NaN], color=c_model["SB"], label=" " * names["SB"], linestyle="-", linewidth=3, zorder=1)
    axs[2,1].legend(loc = "upper right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    for ax in axs[1, :]
        ax.axvline(x_fil[1], linestyle="--", linewidth=2, color="orange")
        ax.axvline(-x_fil[1], linestyle="--", linewidth=2, color="orange")
    end

    for ax in axs[2, :]
        ax.scatter(x_fil[1], x_fil[3], marker="x", s=200, color="orange", linewidth=3)
        ax.scatter(x_fil[1], x_fil[3], marker="o", s=300, edgecolor="orange", facecolor="none", linewidth=3)
    end

    return fig, axs
end

"Make plot showing radius, speed, mass, momentum per SN and kinetic energy efficiencies as a function of time for different models in a structured medium."
function create_figure_filaments(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, model_names::Vector{String}, 
                                 names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, tlim::Vector{<:Real})
    # create figure & axes
    fig = figure(figsize=(18, 10), num=1, clear=true)
    axs = fig.subplots(ncols=3, nrows=2)
    subplots_adjust(wspace=0.25, hspace=0.0)

    # reduce margins
    margins(0,0)

    # plot results
    for model in model_names
        # radius
        axs[1,1].plot(time, models[model][:r][:, 1], color=c_model[model], linestyle="--", linewidth=3, zorder=1)
        axs[1,1].plot(time, models[model][:r_vol][:, 1], color=c_model[model], linewidth=3, zorder=1)
        axs[1,1].fill_between(time, models[model][:r_vol][:, 2], models[model][:r_vol][:, 3], color=c_model[model], linewidth=0, zorder=0, alpha=0.15)

        # velocity
        axs[1,2].plot(time, models[model][:v], color=c_model[model], linewidth=3, zorder=1)
        axs[1,2].plot(time, models[model][:v_t], color=c_model[model], linestyle="--", linewidth=3, zorder=1)

        # Mass
        axs[1,3].plot(time, models[model][:M], color=c_model[model], linewidth=3, zorder=1, label=" " * names[model])

        # Momentum
        axs[2,1].plot(time, models[model][:p], color=c_model[model], linewidth=3, zorder=1)
        axs[2,1].plot(time, models[model][:p_t], color=c_model[model], linewidth=3, zorder=1)

        # Kinetic Energy
        axs[2,2].plot(time, models[model][:f_kin], color=c_model[model], linestyle="-", linewidth=3, zorder=1)
        axs[2,2].plot(time, models[model][:f_kin_t], color=c_model[model], linestyle="--", linewidth=3, zorder=1)
        axs[2,2].plot(time, models[model][:f_pot], color=c_model[model], linestyle="-.", linewidth=3, zorder=1)

        # density
        axs[2,3].plot(time, models[model][:n_H], color=c_model[model], linewidth=3, zorder=1)
    end

    # y-axis labels
    axs[1,1].set_ylabel(L"\mathrm{Shock \  Radius \ [pc]}", fontsize=fontsize)
    axs[1,2].set_ylabel(L"\mathrm{Shock \  Speed \ [km/s]}", fontsize=fontsize)
    axs[1,3].set_ylabel(L"\mathrm{Mass \ [M_{\odot} / sr]}", fontsize=fontsize)
    axs[2,1].set_ylabel(L"\mathrm{Momentum \ per \ SN \ [M_{\odot}/ sr \ km/s]}", fontsize=fontsize)  
    axs[2,2].set_ylabel(L"\mathrm{Energy \ Efficiency}", fontsize=fontsize)
    axs[2,3].set_ylabel(L"\mathrm{Density \ [cm^{-3}]}", fontsize=fontsize)

    # add figure handles
    axs[1,1].text(0.90, 0.05, "(a)", transform=axs[1,1].transAxes, fontsize=fontsize)
    axs[1,2].text(0.90, 0.90, "(b)", transform=axs[1,2].transAxes, fontsize=fontsize)
    axs[1,3].text(0.05, 0.90, "(c)", transform=axs[1,3].transAxes, fontsize=fontsize)
    axs[2,1].text(0.90, 0.05, "(d)", transform=axs[2,1].transAxes, fontsize=fontsize)
    axs[2,2].text(0.05, 0.05, "(e)", transform=axs[2,2].transAxes, fontsize=fontsize)
    axs[2,3].text(0.05, 0.90, "(f)", transform=axs[2,3].transAxes, fontsize=fontsize)

    # decorate axes
    for ax in axs
        # set axis scales
        ax.set_xscale("log");ax.set_yscale("log")
        ax.set_xlim(tlim...)

        # size up tick parameters
        ax.tick_params(axis="both", labelsize=fontsize)
    end

    # remove x-tick labels in top row
    for i in 1:3
        # no x-tick labels for the top row
        axs[1, i].set_xticklabels([])

        # axis label for bottom row
        axs[2,i].set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    end

    # radius panel
    axs[1,1].plot([NaN], [NaN], color="black", label=L" \mathrm{R_{vol}}", linestyle="-", linewidth=1.5, zorder=1)
    axs[1,1].plot([NaN], [NaN], color="black", label=L" \mathrm{R_{ellipsoid}}", linestyle="--", linewidth=1.5, zorder=1)
    axs[1,1].legend(loc = "upper left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # speed panel
    axs[1,2].plot([NaN], [NaN], color="black", label=" Expansion", linestyle="-", linewidth=1.5, zorder=1)
    axs[1,2].plot([NaN], [NaN], color="black", label=" Tangential", linestyle="--", linewidth=1.5, zorder=1)
    axs[1,2].legend(loc = "lower left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # momentum panel
    axs[2,1].plot([NaN], [NaN], color="black", label=" Expansion", linestyle="-", linewidth=1.5, zorder=1)
    axs[2,1].plot([NaN], [NaN], color="black", label=" Tangential", linestyle="--", linewidth=1.5, zorder=1)
    axs[2,1].legend(loc = "upper left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # energy fractions panel
    axs[2,2].plot([NaN], [NaN], color="black", label=" Expansion", linestyle="-", linewidth=1.5, zorder=1)
    axs[2,2].plot([NaN], [NaN], color="black", label=" Tangential", linestyle="--", linewidth=1.5, zorder=1)
    axs[2,2].plot([NaN], [NaN], color="black", label=" Potential", linestyle="-.", linewidth=1.5, zorder=1)
    axs[2,2].legend(loc = "upper right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # legend items
    axs[1,3].plot([NaN], [NaN], color="black", label=" Uniform Medium", linestyle=":", linewidth=3, zorder=1)
    axs[1,3].legend(loc = "lower right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, axs
end

"Plot showing the expansion speed and the mass-weighted expansion speed for a structured medium."
function plot_expansion_speed_filaments(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, 
                                        models_uniform::Dict{String, Dict{Symbol, Vector{Float64}}}, model_names::Vector{String}, 
                                        names::Dict{String, <:AbstractString}, c_model::Dict{String, String})
    # create figure & axes
    fig = figure(figsize=(8, 8), num=1, clear=true)
    ax = gca()

    # reduce margins
    margins(0,0)

    for model in model_names
        # velocity from momentum
        ax.plot(time, models[model][:v], color=c_model[model], linestyle="-", linewidth=3, zorder=1)
        ax.plot(time, models_uniform[model][:v], color=c_model[model], linestyle=":", linewidth=3, zorder=1)

        # get gradient from rate of change of volume
        t_grad, v_grad = gradient(time, models[model][:r_vol][:, 1])
        ax.plot(t_grad, v_grad * pc / Myr / km_s, color=c_model[model], linestyle="--", linewidth=3, zorder=1)
    end

    # decorate axis
    # set axis scales
    ax.set_xscale("log");ax.set_yscale("log")

    # Axis labels
    ax.set_ylabel(L"\mathrm{Shock \  Speed \ [km/s]}", fontsize=fontsize)
    ax.set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)

    # legend panel
    ax.plot([NaN], [NaN], color="black", label=L"\, \mathrm{v_{M}}", linestyle="-", linewidth=3, zorder=1)
    ax.plot([NaN], [NaN], color="black", label=L"\,  \mathrm{v_{Vol}}", linestyle="--", linewidth=3, zorder=1)
    ax.plot([NaN], [NaN], color=c_model["SN"], label=" " * names["SN"], linestyle="-", linewidth=3, zorder=1)
    ax.plot([NaN], [NaN], color=c_model["SB"], label=" " * names["SB"], linestyle="-", linewidth=3, zorder=1)
    ax.plot([NaN], [NaN], color="black", label=" Uniform Medium", linestyle=":", linewidth=3, zorder=1)
    ax.legend(loc = "lower left", frameon=false, handlelength=1.5, fontsize=0.8fontsize)

    return fig, ax
end

"Plot showing geometry tracks for different models in a structured medium."
function plot_geometry_tracks_filaments(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, model_names::Vector{String}, 
                                names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, times_orb::Union{Vector{<:Real}, AbstractRange})
    # get orbital times
    i_orb = get_timepoints(time, times_orb)
    t_max = times_orb[end]
    i_t   = time .< t_max

    # create figure & axes
    fig = figure(figsize=(8, 8), num=1, clear=true)

    # create geometry space
    ax = gca()
    create_geometry_space!(ax)

    # place colorbars
    norm_t = matplotlib.colors.Normalize(vmin=times_orb[1], vmax=times_orb[end], clip=true)
    colors_time = ColorSchemes.managua[LinRange(0, 1.0, 256)]
    cmap_time = ColorScheme(colors_time)
    sm   = plt.cm.ScalarMappable(cmap=ColorMap(colors_time), norm=norm_t)
    s    = ax.get_position()
    cbaxes = fig.add_axes([s.x0 + 0.03 * s.width, s.y0 + 0.3 * s.height, 0.03 * s.width, 0.4 * s.height])
    cb = colorbar(sm, cax=cbaxes, orientation="vertical")
    cb.set_label(label=L"\mathrm{Time \ [Myr]}", size=0.7fontsize)
    cb.ax.tick_params(axis="y", labelsize=0.7fontsize)

    # reduce margins
    margins(0,0)

    # get time-> color mapping
    c_time = [get_color(cmap_time, t, norm=norm_t) for t in times_orb]

    # plot model
    for model in model_names
        ax.plot(models[model][:major_ratio][i_t],   models[model][:minor_ratio][i_t],   color=c_model[model], linewidth=3, zorder=1, label=" " * names[model])
        ax.scatter(models[model][:major_ratio][i_orb], models[model][:minor_ratio][i_orb], c=c_time, edgecolors=c_model[model], marker=".", s=75, zorder=2)
    end

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[end-1:end], labels[end-1:end], loc = "upper left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, ax
end