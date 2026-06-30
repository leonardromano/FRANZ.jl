############################################################################################################################################
# Functions for axis geometry and decoration

"Return a MultipleLocator object with less verbose syntax"
MultipleLocator(number::Number) = matplotlib.ticker.MultipleLocator(number)

"Symmetrically set the major and minor tick locators of both axes for a given axis object."
function specify_ticks(axis; major::Float64=0.5, minor::Float64=0.1)
    axis.xaxis.set_major_locator(MultipleLocator(major)) # Major ticks: multiples of major kpc
    axis.xaxis.set_minor_locator(MultipleLocator(minor)) # Major ticks: multiples of minor kpc
    axis.yaxis.set_major_locator(MultipleLocator(major)) # Major ticks: multiples of major kpc
    axis.yaxis.set_minor_locator(MultipleLocator(minor)) # Major ticks: multiples of minor kpc
end

"compute a reasonable scale for major and minor ticks"
function get_tick_scale(L::Real)
    mag = floor(Int, log10(L))
    s   = 10.0^mag

    if L / s < 5
        return s, 0.2 * s, 0.5 * s / L
    else
        return 5 * s, s, s / L
    end
end

"Return the number truncated to the level of significant figures"
function significant_figures(number::Float64; digits::Integer=2)
    OoM = (number == 0.0) ? 0 : floor(Int, log10(abs(number)))
    return (OoM > 0) ? round(number, digits=digits) : round(number, digits=digits-OoM)
end

"Return a color given a colorscheme and a value"
function get_color(colorscheme, value; norm=nothing)
    norm = isnothing(norm) ? f(x) = x : norm 
    color = colorscheme[norm(value)]
    return (red(color), green(color), blue(color))
end

"function for adding timescale annotations to the plots"
function add_timescale_annotation(ax, t::Real; ls_ann="-.", lw_ann=2, label=nothing, ylabel::Real=0.0, color="black")
    if !isnothing(label)
        line = ax.plot([t-SMALL, t+SMALL], ylabel * [10^0.1, 10^-0.1], color=color, linestyle=ls_ann, zorder=-1, linewidth=lw_ann, label=label)
        ll.labelLines(line, align=true, xvals=[t], fontsize=0.7fontsize)
    end

    # plot vertical line
    ax.axvline(t, color=color, linestyle=ls_ann, zorder=-1, linewidth=lw_ann)
end

"function for adding timescale annotations to the plots"
function add_horizontal_annotation(ax, y::Real; ls_ann="-.", lw_ann=2, label=nothing, xlabel::Real=0.0, color="black")
    if !isnothing(label)
        line = ax.plot(xlabel * [10^0.1, 10^-0.1], [y-SMALL, y+SMALL], color=color, linestyle=ls_ann, zorder=-1, linewidth=lw_ann, label=label)
        ll.labelLines(line, align=true, xvals=[xlabel], fontsize=0.7fontsize)
    end

    # plot horizontal line
    ax.axhline(y, color=color, linestyle=ls_ann, zorder=-1, linewidth=lw_ann)
end

"function for adding line annotations to the plots"
function add_line_annotation(ax, x::Vector{<:Real}, y::Vector{<:Real}; ls_ann="-.", lw_ann=2, label=nothing, ilabel=2, color="black")
    if isnothing(label)
        ax.plot(x, y, color=color, linestyle=ls_ann, zorder=2, linewidth=lw_ann)
    else
        ilabels = max(1, ilabel):min(length(x), ilabel+1)

        line = ax.plot((1-SMALL) * x[ilabels], y[ilabels], color=color, linestyle=ls_ann, zorder=2, linewidth=0.1lw_ann, label=label)
        ll.labelLines(line, align=true, xvals=[x[ilabel]], fontsize=0.7fontsize)
        ax.plot(x[1:ilabel], y[1:ilabel], color=color, linestyle=ls_ann, zorder=2, linewidth=lw_ann)
    end
end

############################################################################################################################################
# Function for setting the style defaults

"Set the default plotting parameters to ensure a uniform style"
function set_plot_style()
    # linewidth
    matplotlib.rc("axes", linewidth=1)
    matplotlib.rc("lines", linewidth=1)

    ## ticks
    # x-ticks
    matplotlib.rc("xtick", top=true, bottom=true, direction="in")
    matplotlib.rc("xtick.major", size=6, width=1, pad=4)
    matplotlib.rc("xtick.minor", size=3, width=1, pad=4, visible=true)

    # y-ticks (the same)
    matplotlib.rc("ytick", left=true, right=true, direction="in")
    matplotlib.rc("ytick.major", size=6, width=1, pad=4)
    matplotlib.rc("ytick.minor", size=3, width=1, pad=4, visible=true)
    
    # legend
    matplotlib.rc("legend", frameon=false, handletextpad=0.0)

    # figure
    matplotlib.rc("figure", figsize=(4, 4))
    matplotlib.rc("savefig", dpi=300, format="pdf")
end

############################################################################################################################################
# Functions for setting up layouts for specific plots

"Create the backbone of a geometry space plot."
function create_geometry_space!(ax)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)

    # plot special regions (differing from Colman+24; van der Wel+14)
    ax.fill_between([0, 1], [0, 1], [1,1], color="red", alpha=0.3, zorder=1) # forbidden region
    ax.plot([0,1], [0,1], color="k", linewidth=2, label="a=b", zorder=0)
    ax.plot([0.6667, 1.0], [0.6667, 0.6667], color="k", linewidth=2, linestyle=":", zorder=-1, label="")
    ax.plot([0.6667, 0.6667], [0, 0.6667], color="k", linewidth=2, linestyle=":", zorder=-1, label="")
    ax.plot([0.6667, 1], [0.44444, 0.6667], color="k", linewidth=2, linestyle=":", label="a=2/3 b", zorder=-1)
    ll.labelLines(ax.get_lines(), align=true, xvals=[0.5, 5/6], zorder=0, fontsize=fontsize)

    # annotate each region
    ax.annotate("P", xy=(0.3333, 0.166667), fontsize=fontsize, weight="bold") # Prolate
    ax.annotate("S", xy=(0.8333, 0.75), fontsize=fontsize, weight="bold")     # Spheroidal
    ax.annotate("O", xy=(0.8333, 0.3333), fontsize=fontsize, weight="bold")   # Oblate
    ax.annotate("OS", xy=(0.75, 0.6), fontsize=fontsize, weight="bold")       # Oblate Spheroidal

    ax.set_xlabel("Semi-Major-to-Major-Ratio b/c", fontsize=fontsize)
    ax.set_ylabel("Minor-to-Major-Ratio a/c", fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize)
end

