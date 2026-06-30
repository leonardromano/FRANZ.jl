############################################################################################################################################################
# functions to run a model and return its output in a nice format for plotting

"Return the output of the model for a 3D medium"
function run_model_shear(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0, potential::Function=z->0.0, Nside::Integer=4)
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, map=HealpixMap{Float64, NestedOrder}(Nside))

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 * input2func(model.E_51)

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

        # get galactocentric coordinates
        z = [xi[3] for xi in x[it, :]]
        R = [norm(xi[1:2]) for xi in x[it, :]]

        # get SNR radius
        r_shear[it, 1] = prod(abc)^(1/3)
        r_shear[it, 2] = minimum(abc)
        r_shear[it, 3] = maximum(abc)

        minor_ratio[it] = abc[3] / abc[1]
        major_ratio[it] = abc[2] / abc[1]

        # get angles
        e_φ = normalize(cross([0, 0, 1], x_s[it, :]))
        e_R = normalize(cross([0, 0, 1], e_φ))

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
        x_c = env_shear.x_c(t_shear[it])
        Mass_shear[it]   = sum(M[it, :]) * dΩ
        Volume_shear[it] = sum([dot(x[it, ihp] - x_c, n[it, ihp]) / 3 for ihp in eachindex(x[it, :])]) * dΩ
        Area_shear[it]   = sum([norm(dA) for dA in n[it, :]]) * dΩ
        r_vol[it, 1]     = (3 * Volume_shear[it] / 4π)^(1/3)
        r_c              = [norm(xi - x_c) for xi in x[it, :]]
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
        p_shear[it] = p_rad / Num_SN(t[it, 1])
        p_t[it]     = p_perp / Num_SN(t[it, 1])

        # get kinetic energy
        f_kin_rad_shear[it] = 0.5 * sum(M[it, :] .* v_rad.^2) * dΩ * km_s^2 * Msol / E_SN(t[it, 1])
        f_kin_t_shear[it]   = 0.5 * sum(M[it, :] .* v_t.^2) * dΩ * km_s^2 * Msol / E_SN(t[it, 1])

        # get potential energy
        f_pot_shear[it] = sum(M[it, :] .* potential.(R)) * Msol * dΩ / E_SN(t[it, 1])
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
    output[:A]       = Area_shear
    output[:p]       = p_shear
    output[:p_t]     = p_t
    output[:f_kin]   = f_kin_rad_shear
    output[:f_kin_t] = f_kin_t_shear
    output[:f_pot]   = f_pot_shear

    return output
end

"Return the output of the model for a slice through a 3D medium"
function run_model_shear_slice(model::Model, environment::Environment, time::Union{AbstractRange{<:Real}, Vector{<:Real}}; N_SN=1.0,
    φ_range::Union{AbstractRange{<:Real}, Vector{<:Real}}=LinRange(0, 2π, 5))
    # run model
    t, x, v, n, M, f = numerical_solution(time, model=model, environment=environment, cosθ=0.0, ϕ=φ_range)

    # number of SNe as a function of time
    Num_SN = input2func(N_SN)
    E_SN   = 1e51 / 4π * input2func(model.E_51)

    output = Dict{Symbol, Any}()
    output[:t] = t
    output[:x] = x
    output[:n] = n
    output[:v] = v
    output[:p] = M .* v ./ Num_SN.(t)

    return output
end

############################################################################################################################################################
# functions for plotting the output

"Make plot showing radius, speed, mass, momentum per SN and kinetic energy efficiencies as a function of time for different models in a shearing medium."
function create_figure_shear(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, model_names::Vector{String}, 
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
        axs[1,3].plot(time, models[model][:M], color=c_model[model], linewidth=3, zorder=1)

        # Momentum
        axs[2,1].plot(time, models[model][:p], color=c_model[model], linewidth=3, zorder=1)
        axs[2,1].plot(time, models[model][:p_t], color=c_model[model], linewidth=3, zorder=1)

        # Kinetic Energy
        axs[2,2].plot(time, models[model][:f_kin], color=c_model[model], linestyle="-", linewidth=3, zorder=1)
        axs[2,2].plot(time, models[model][:f_kin_t], color=c_model[model], linestyle="--", linewidth=3, zorder=1)
        axs[2,2].plot(time, models[model][:f_kin] + models[model][:f_kin_t], color=c_model[model], linestyle="-.", linewidth=3, zorder=1)

        # Legend panel
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
    axs[1,1].text(0.05, 0.05, "(a)", transform=axs[1,1].transAxes, fontsize=fontsize)
    axs[1,2].text(0.90, 0.90, "(b)", transform=axs[1,2].transAxes, fontsize=fontsize)
    axs[1,3].text(0.05, 0.90, "(c)", transform=axs[1,3].transAxes, fontsize=fontsize)
    axs[2,1].text(0.90, 0.05, "(d)", transform=axs[2,1].transAxes, fontsize=fontsize)
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
    axs[2,2].plot([NaN], [NaN], color="black", label=" Total Kin.", linestyle="-.", linewidth=1.5, zorder=1)
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

"Plot showing geometry tracks for different models in a shearing medium."
function plot_geometry_tracks_shear(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, model_names::Vector{String}, 
                                names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, times_orb::Union{Vector{<:Real}, AbstractRange}, t_orb::Real)
    # get orbital times
    i_orb = get_timepoints(time, times_orb * t_orb)
    i_def = Dict("SN" => findmin(abs.(t_orb/8π .- time))[2], "SB" => findmin(abs.(t_orb/4π .- time))[2], "StB" => findmin(abs.(t_orb/4π .- time))[2])

    # create figure & axes
    fig = figure(figsize=(8, 8), num=1, clear=true)

    # create geometry space
    ax = gca()
    create_geometry_space!(ax)

    # place colorbars
    norm_t = matplotlib.colors.Normalize(vmin=times_orb[1], vmax=times_orb[end], clip=true)
    colors_time = ColorSchemes.managua[LinRange(0.0, 1.0, 256)]
    cmap_time = ColorScheme(colors_time)
    sm   = plt.cm.ScalarMappable(cmap=ColorMap(colors_time), norm=norm_t)
    s    = ax.get_position()
    cbaxes = fig.add_axes([s.x0 + 0.03 * s.width, s.y0 + 0.3 * s.height, 0.03 * s.width, 0.4 * s.height])
    cb = colorbar(sm, cax=cbaxes, orientation="vertical")
    cb.set_label(label=L"\mathrm{Time \ [t_{orb}]}", size=0.7fontsize)
    cb.ax.tick_params(axis="y", labelsize=0.7fontsize)
    cb.ax.yaxis.set_ticks([times_orb...], minor=false)

    # reduce margins
    margins(0,0)

    # get time-> color mapping
    c_time = [get_color(cmap_time, t, norm=norm_t) for t in times_orb]

    # plot model
    for model in model_names
        ax.plot(models[model][:major_ratio][1:i_orb[5]],   models[model][:minor_ratio][1:i_orb[5]],   color=c_model[model], linewidth=3, zorder=1, label=" " * names[model])
        ax.plot(models[model][:major_ratio][i_orb[5]:end], models[model][:minor_ratio][i_orb[5]:end], color=c_model[model], linewidth=3, zorder=1, linestyle=":", alpha=0.25)
        ax.scatter(models[model][:major_ratio][i_orb], models[model][:minor_ratio][i_orb], c=c_time, edgecolors=c_model[model], marker=".", s=75, zorder=2)
        ax.scatter(models[model][:major_ratio][i_def[model]], models[model][:minor_ratio][i_def[model]], color="gold", edgecolors=c_model[model], marker="*", s=100, zorder=2)
    end

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[end-2:end], labels[end-2:end], loc = "upper left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, ax
end

"Plot showing the pitch angle as a function for different models in a shearing medium."
function plot_pitch_angles_shear(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Array{Float64}}}, model_names::Vector{String}, 
                                names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, times_orb::Union{Vector{<:Real}, AbstractRange}, t_orb::Real)
    # get orbital times
    i_orb = get_timepoints(time, times_orb * t_orb)
    i_def = Dict("SN" => findmin(abs.(t_orb/8π .- time))[2], "SB" => findmin(abs.(t_orb/4π .- time))[2], "StB" => findmin(abs.(t_orb/4π .- time))[2])

    # create figure & axes
    fig = figure(figsize=(8, 8), num=1, clear=true)

    # create axis
    ax = gca()
    ax.set_xlim(0.25 * t_orb / (8π), t_shear[end])
    ax.set_xscale("log")
    ax.set_ylim(0, 50)

    ax.set_ylabel("Pitch Angle [deg]", fontsize=fontsize)
    ax.set_xlabel(L"\mathrm{Time \ [Myr]}", fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize)

    # reduce margins
    margins(0,0)

    first_half = falses(size(time))
    first_half[1:i_orb[5]] .= true

    # plot model
    for model in model_names
        # only show points where a/c < 2/3
        deformed   = models[model][:minor_ratio].< 2/3
        
        # major axis
        ax.plot(t_shear[first_half .& deformed], models[model][:ϕ_major][first_half .& deformed], color=c_model[model], label=" " * names[model], linewidth=3, zorder=1, linestyle="-")
        ax.plot(t_shear[first_half .& .!(deformed)], models[model][:ϕ_major][first_half .& .!(deformed)], color=c_model[model], linewidth=3, zorder=1, linestyle="--", alpha=0.25)
        ax.plot(t_shear[i_orb[5]:end], models[model][:ϕ_major][i_orb[5]:end], color=c_model[model], linewidth=3, zorder=1, linestyle=":", alpha=0.25)
        ax.scatter(t_shear[i_orb], models[model][:ϕ_major][i_orb], color=c_model[model], edgecolors=c_model[model], marker=".", s=75, zorder=2)
        ax.scatter(t_shear[i_def[model]], models[model][:ϕ_major][i_def[model]], color="gold", edgecolors=c_model[model], marker="*", s=100, zorder=2)
    end

    # legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[end-2:end], labels[end-2:end], loc = "upper right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, ax
end

"Plot showing the contour lines of the shock surface at different times for different models. Also show velocity vectors as arrows."
function plot_contours_shear(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Any}}, environment::Environment, model_names::Vector{String}, 
                             names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, times_frame::Vector{<:Real}, t_orb::Real; 
                             L_frame::Vector{<:Vector{<:Real}}=[[125, 1125, 2000], [625, 1e4, 1.75e4]], 
                             vel_scale::Vector{<:Vector{<:Real}}=[[500, 15, 15], [4000, 200, 150]],
                             N_arrow_sampling::Integer=4)

    # create figure
    fig = figure(figsize=(18, 12), num=1, clear=true)
    axs = fig.subplots(ncols=3, nrows=2)
    subplots_adjust(wspace=0.01, hspace=0.01)

    # for which subset to plot arrows?
    Nφ        = size(models["StB"][:x])[2]
    Δi_sample = max(div(Nφ, N_arrow_sampling), 1)
    samples   = 1:Δi_sample:Nφ

    # get orbital times for each panel
    i_orb = get_timepoints(time, times_frame * t_orb)

    for i_time in eachindex(times_frame)
        # get axis and center
        i_out = i_orb[i_time]
        x0 = env_shear.x_c(time[i_out])
        e_rot = normalize(env_shear.v_ext(x0, time[i_out])[1:2])
        e_R   = normalize(x0[1:2])

        # top row: SN & SB
        ax = axs[1, i_time]

        # plot explosion center
        ax.scatter(x0[1], x0[2], color="orange", marker="*", s=150)
        
        for model in ["SN", "SB"]
            ax.plot([pos[1] for pos in models[model][:x][i_out, :]], [pos[2] for pos in models[model][:x][i_out, :]], color=c_model[model], linewidth=3)
            ax.quiver([pos[1] for pos in models[model][:x][i_out, samples]], [pos[2] for pos in models[model][:x][i_out, samples]], 
            [pos[1] for pos in models[model][:v][i_out, samples]], [pos[2] for pos in models[model][:v][i_out, samples]], width=0.005, color=c_model[model], scale=vel_scale[1][i_time])
        end

        ax.set_xlim(x0[1]-0.5L_frame[1][i_time], x0[1]+0.5L_frame[1][i_time])
        ax.set_ylim(x0[2]-0.5L_frame[1][i_time], x0[2]+0.5L_frame[1][i_time])

        # decorate axes
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.annotate(L"\mathrm{t = %$(times_frame[i_time]) \ t_{orb}}", xy=x0[1:2] + [-0.45, 0.425] * L_frame[1][i_time], fontsize=fontsize, color="black")

        x_ruler = x0[1] .+ [-0.45, -0.45 + 0.2] * L_frame[1][i_time]; y_ruler = x0[2] .+ [-0.45, -0.45] * L_frame[1][i_time]
        ax.plot(x_ruler, y_ruler, color="black", linewidth=2)

        shift = i_time == 1 ? 0.03 : 0.0125
        ax.annotate("$(round(Int, 0.2 * L_frame[1][i_time])) pc", xy=(x_ruler[1]+ shift * L_frame[1][i_time], mean(y_ruler)+0.025 * L_frame[1][i_time]), fontsize=fontsize, color="black")  

        # place compass
        # rotation direction
        dx_arrow = 0.1 * L_frame[1][i_time] * e_rot
        x_arrow  = x0[1:2] + [0.4, 0.4] * L_frame[1][i_time] - 0.5 * dx_arrow
        arrow_width = 0.005 * L_frame[1][i_time]
        ax.arrow(x_arrow[1], x_arrow[2], dx_arrow[1], dx_arrow[2], color="crimson", width=arrow_width, head_width=5 * arrow_width)

        # galactic center
        dx_arrow = -0.1 * L_frame[1][i_time] * e_R
        x_arrow  = x0[1:2] + [0.4, 0.4] * L_frame[1][i_time] - 0.5 * dx_arrow
        ax.arrow(x_arrow[1], x_arrow[2], dx_arrow[1], dx_arrow[2], color="cornflowerblue", width=arrow_width, head_width=5 * arrow_width)
        
        # plot StB separately
        ax = axs[2, i_time]

        # plot explosion center
        ax.scatter(x0[1], x0[2], color="orange", marker="*", s=150)

        ax.plot([pos[1] for pos in models["StB"][:x][i_out, :]], [pos[2] for pos in models["StB"][:x][i_out, :]], color=c_model["StB"], linewidth=3)
        ax.quiver([pos[1] for pos in models["StB"][:x][i_out, samples]], [pos[2] for pos in models["StB"][:x][i_out, samples]], 
            [pos[1] for pos in models["StB"][:v][i_out, samples]], [pos[2] for pos in models["StB"][:v][i_out, samples]], width=0.005, color=c_model["StB"], scale=vel_scale[2][i_time])

        # set axis limits
        ax.set_xlim(x0[1]-0.5L_frame[2][i_time], x0[1]+0.5L_frame[2][i_time])
        ax.set_ylim(x0[2]-0.5L_frame[2][i_time], x0[2]+0.5L_frame[2][i_time])

        # decorate axes
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        ax.annotate(L"\mathrm{t = %$(times_frame[i_time]) \ t_{orb}}", xy=x0[1:2] + [-0.45, 0.425] * L_frame[2][i_time], fontsize=fontsize, color="black")

        x_ruler = x0[1] .+ [-0.45, -0.45 + 0.2] * L_frame[2][i_time]; y_ruler = x0[2] .+ [-0.45, -0.45] * L_frame[2][i_time]
        ax.plot(x_ruler, y_ruler, color="black", linewidth=2)

        shift = i_time == 1 ? 0.0125 : 0.0
        ax.annotate("$(round(Int, 0.2 * L_frame[2][i_time])) pc", xy=(x_ruler[1] + shift * L_frame[2][i_time], mean(y_ruler)+0.025 * L_frame[2][i_time]), fontsize=fontsize, color="black")

        # place compass
        # rotation direction
        dx_arrow = 0.1 * L_frame[2][i_time] * e_rot
        x_arrow  = x0[1:2] + [0.4, 0.4] * L_frame[2][i_time] - 0.5 * dx_arrow
        arrow_width = 0.005 * L_frame[2][i_time]
        ax.arrow(x_arrow[1], x_arrow[2], dx_arrow[1], dx_arrow[2], color="crimson", width=arrow_width, head_width=5 * arrow_width)

        # galactic center
        dx_arrow = -0.1 * L_frame[2][i_time] * e_R
        x_arrow  = x0[1:2] + [0.4, 0.4] * L_frame[2][i_time] - 0.5 * dx_arrow
        ax.arrow(x_arrow[1], x_arrow[2], dx_arrow[1], dx_arrow[2], color="cornflowerblue", width=arrow_width, head_width=5 * arrow_width)
    end

    # add legend
    axs[1,1].plot([NaN], [NaN], color=c_model["SN"], label=" " * names["SN"], linestyle="-", linewidth=3, zorder=1)
    axs[1,1].legend(loc = "lower right", frameon=false, handlelength=1, fontsize=0.7fontsize)
    axs[1,2].plot([NaN], [NaN], color=c_model["SB"], label=" " * names["SB"], linestyle="-", linewidth=3, zorder=1)
    axs[1,2].legend(loc = "lower right", frameon=false, handlelength=1, fontsize=0.7fontsize)
    axs[2,1].plot([NaN], [NaN], color=c_model["StB"], label=" " * names["StB"], linestyle="-", linewidth=3, zorder=1)
    axs[2,1].legend(loc = "lower right", frameon=false, handlelength=1, fontsize=0.7fontsize)

    return fig, axs
end

"Plot showing the different momentum components as a function of time for different models in a shearing medium to highlight the effect of epicycles."
function plot_momentum_decomposition(time::Union{Vector{<:Real}, AbstractRange}, models::Dict{String, Dict{Symbol, Any}}, models_uniform::Dict{String, Dict{Symbol, Vector{Float64}}}, 
                                     environment::Environment, model_names::Vector{String}, names::Dict{String, <:AbstractString}, c_model::Dict{String, String}, 
                                     t_end::Real, t_orb::Real, φ_range::Union{Vector{<:Real}, AbstractRange})
    # Create single figure panel
    fig = figure(figsize=(18, 6), num=1, clear=true)
    axs = fig.subplots(ncols=3, nrows=1)
    subplots_adjust(wspace=0.01, hspace=0.00)

    # reduce margins
    margins(0,0)

    # get orbital times
    samples = 0.0:(π/2):(π/2)
    i_φs = get_timepoints(φ_range, samples)

    # coloring of angles
    norm_φ = matplotlib.colors.Normalize(vmin=0.0, vmax=2π, clip=true)
    colors_φ = ColorSchemes.phase[LinRange(0.0, 1.0, 256)]
    cmap_φ = ColorScheme(colors_φ)

    # decorate the axes
    for ax in axs
        ax.set_xlim(0.01, t_end)
        ax.set_ylim(-1.45, 1.45)
        ax.set_xlabel(L"\mathrm{Time \ [t_{\kappa}]}", fontsize=fontsize)
        ax.tick_params(axis="both", labelsize=fontsize)

        ax.axhline(0.0, linestyle=":", color="gray")
        ax.axhline(1.0, linestyle="-", color="gray")
        ax.axhline(-1.0, linestyle="-", color="gray")
        ax.axhline(1 / sqrt(2), linestyle="--", color="gray")
        ax.axhline(-1 / sqrt(2), linestyle="--", color="gray")
    end

    # plot model lines
    for (i, model) in enumerate(model_names)
        # get axis
        ax = axs[i]
        ax.set_title(names[model], fontsize=fontsize)
        
        for i_φ in i_φs
            # get velocity decomposition
            x = models[model][:x][:, i_φ]
            v = models[model][:p][:, i_φ] ./ models_uniform[model][:p]

            e_R = [normalize([pos[1:2]..., 0.0]) for pos in x]
            e_φ = [cross([0, 0, 1], e) for e in e_R]

            v_R = [e'vel for (vel, e) in zip(v, e_R)]
            v_φ = [e'vel for (vel, e) in zip(v, e_φ)]
            v_norm = [norm(vel) for vel in v]

            # plot velocities
            c_φ = get_color(cmap_φ, φ_range[i_φ], norm=norm_φ)
            ax.plot(sqrt(2) * time / t_orb, v_norm, color=c_φ, linewidth=3, zorder=1, linestyle="-", alpha=0.5)
            ax.plot(sqrt(2) * time / t_orb, v_R, color=c_φ, linewidth=3, zorder=1, linestyle="--", alpha=0.5)
            ax.plot(sqrt(2) * time / t_orb, v_φ, color=c_φ, linewidth=3, zorder=1, linestyle=":", alpha=0.5)
        end
    end

    # add some stuff to the first axis
    ax = axs[1]
    ax.set_ylabel(L"\mathrm{Momentum \ [p_{inj}(t)]}", fontsize=fontsize)

    # plot model
    for i_φ in i_φs
        # model
        phase = LinRange(0.0, t_end, 1000)

        v_ampl = 1 / sqrt(1 + cos(φ_range[i_φ])^2)
        v_R_model = sqrt(2) * v_ampl * cos.(2π * phase.- φ_range[i_φ])
        v_φ_model = - v_ampl * sin.(2π * phase .- φ_range[i_φ])
        v_model   = sqrt.(v_R_model.^2 + v_φ_model.^2)
        
        c_φ = get_color(cmap_φ, φ_range[i_φ], norm=norm_φ)
        ax.plot(phase, v_R_model, color=c_φ, linestyle="--", alpha=1, linewidth=1.5)
        ax.plot(phase, v_φ_model, color=c_φ, linestyle=":", alpha=1, 1.5)
        ax.plot(phase, v_model, color=c_φ, linestyle="-", alpha=1, 1.5)
    end
    
    # legend labels
    ax.plot([NaN], [NaN], color="gray", label=" |p|", linewidth=3, zorder=1, linestyle="-")
    ax.plot([NaN], [NaN], color="gray", label=L" \mathrm{p_{R}}", linewidth=3, zorder=1, linestyle="--")
    ax.plot([NaN], [NaN], color="gray", label=L" \mathrm{p_{\varphi}}", linewidth=3, zorder=1, linestyle=":")
    ax.plot([NaN], [NaN], color="black", linestyle="-", alpha=1, linewidth=1.5, label=" Eqs. (2.41) & (2.42)")
    ax.legend(loc = "lower left", frameon=false, handlelength=1, fontsize=0.7fontsize)

    # add some stuff to the second axis
    for i_φ in i_φs
        c_φ = get_color(cmap_φ, φ_range[i_φ], norm=norm_φ)
        axs[2].plot([NaN], [NaN], color=c_φ, label=" φ = $(φ_range[i_φ] / π) π", linewidth=3, zorder=1, linestyle="-", alpha=0.5)
    end
    axs[2].legend(loc = "lower left", frameon=false, handlelength=1, fontsize=0.7fontsize)
    
    # 2nd & 3rd axis: no tick labels
    for i in 2:3
        axs[i].set_yticklabels([])
    end

    return fig, axs
end