using Plots
using DelimitedFiles

# Helper: safely read first two numeric columns from a file
function read_xy_columns(filename::String)
    # Determine delimiter. Default to whitespace/tab, use comma for csv
    delim = endswith(lowercase(filename), ".csv") ? ',' : ' '
    
    # Read file. 'readdlm' handles detecting dimensions
    raw_data = readdlm(filename, delim)
    
    # Heuristic to skip header: check if 1st element (1,1) is a number.
    # If not, assume row 1 is a header.
    start_row = 1
    val_1_1 = raw_data[1, 1]
    
    is_number = (val_1_1 isa Number) || (val_1_1 isa AbstractString && !isnothing(tryparse(Float64, val_1_1)))
    
    if !is_number
        start_row = 2
    end
    
    # Extract
    x_raw = raw_data[start_row:end, 1]
    y_raw = raw_data[start_row:end, 2]
    
    # Convert to floats if needed
    x = Float64.(x_raw)
    y = Float64.(y_raw)
    
    return x, y
end

# Helper: safely read columns 1, 2, and 3 from a file (X, Y1, Y2)
function read_xyy_columns(filename::String)
    # Determine delimiter
    delim = endswith(lowercase(filename), ".csv") ? ',' : ' '
    
    # Read file
    raw_data = readdlm(filename, delim)
    
    # Heuristic to skip header
    start_row = 1
    val_1_1 = raw_data[1, 1]
    is_number = (val_1_1 isa Number) || (val_1_1 isa AbstractString && !isnothing(tryparse(Float64, val_1_1)))
    
    if !is_number
        start_row = 2
    end
    
    # Extract
    x_raw = raw_data[start_row:end, 1]
    y1_raw = raw_data[start_row:end, 2]
    y2_raw = raw_data[start_row:end, 3]
    
    # Convert to floats
    x = Float64.(x_raw)
    y1 = Float64.(y1_raw)
    y2 = Float64.(y2_raw)
    
    return x, y1, y2
end

# Format x tick labels as integers for cleaner axis display
int_tick_label(x) = string(round(Int, x))

"""
    plot_file(filename, title, xlabel, ylabel;
              output_file=nothing, xscale=1.0,
              ticks_from_data=nothing, tick_font::Union{Nothing,Real}=nothing,
              legend_position::Symbol = :none,
              rotate_x_labels::Bool=false)

Reads first two columns of `filename` (x, y) and plots y vs x.

- `xscale`: multiply the x values by this factor before plotting.
- `ticks_from_data`: when `true`, force x ticks to appear at every scaled
  data point (useful when the automatic ticks skip values after scaling).
- `legend_position`: legend placement symbol (default `:none`).
- `rotate_x_labels`: when `true`, rotate x-axis tick labels by 90°.
"""
function plot_file(filename::String, x_label::String, y_label::String;
                   output_file=nothing, xscale::Real=1.0,
                   ticks_from_data::Union{Bool,Nothing}=nothing,
                   tick_font::Union{Nothing,Real}=nothing,
                                     legend_position::Symbol = :none,
                                     rotate_x_labels::Bool=false)
    x, y = read_xy_columns(filename)
    x_scaled = xscale == 1 ? x : x .* xscale
    
    p = plot(x_scaled, y,
        xlabel = x_label,
        ylabel = y_label,
        lw = 2,           # line width
        marker = :circle, # simple marker
        legend = legend_position,
        xrotation = rotate_x_labels ? 90 : 0,
        xformatter = int_tick_label
    )

    use_data_ticks = isnothing(ticks_from_data) ? (xscale != 1) : ticks_from_data
    if use_data_ticks
        ticks = sort(unique(x_scaled))
        labels = int_tick_label.(ticks)
        plot!(p; xticks = (ticks, labels), xlims = (minimum(ticks), maximum(ticks)))
    end

    if !isnothing(tick_font)
        plot!(p; tickfont = font(tick_font))
    end
    
    # Determine output file name
    pdf_file = if isnothing(output_file)
        # Default: use input filename with .pdf extension
        base = splitext(filename)[1]
        "$(base).pdf"
    else
        # Ensure .pdf extension
        endswith(output_file, ".pdf") ? output_file : "$(output_file).pdf"
    end
    
    # Ensure directory exists and save
    dir = dirname(pdf_file)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    savefig(p, pdf_file)
    println("Saved plot to: $pdf_file")
    return p
end

"""
    plot_comparison(file1, file2, legend1, legend2, xlabel, ylabel;
                    output_file=nothing, xscale=1.0,
                    ticks_from_data=nothing, tick_font::Union{Nothing,Real}=nothing,
                    legend_position::Symbol = :best,
                    rotate_x_labels::Bool=false)

Reads two files and plots them on the same figure for comparison.

- `xscale`: multiply every x value by this factor before plotting.
- `ticks_from_data`: when `true`, force x ticks to include all scaled data
  points from both files so the axis labels stay aligned with the data.
- `legend_position`: legend placement symbol (default `:best`).
- `rotate_x_labels`: when `true`, rotate x-axis tick labels by 90°.
"""
function plot_comparison(file1::String, file2::String, legend1::String, legend2::String,
                         x_label::String, y_label::String;
                         output_file=nothing, xscale::Real=2.0,
                         ticks_from_data::Union{Bool,Nothing}=nothing,
                         tick_font::Union{Nothing,Real}=nothing,
                                                 legend_position::Symbol = :best,
                                                 rotate_x_labels::Bool=false)
    x1, y1 = read_xy_columns(file1)
    x2, y2 = read_xy_columns(file2)
    x1_scaled = xscale == 1 ? x1 : x1 .* xscale
    x2_scaled = xscale == 1 ? x2 : x2 .* xscale

    p = plot(
        xlabel = x_label,
        ylabel = y_label,
        legend = legend_position,
        xrotation = rotate_x_labels ? 90 : 0,
        xformatter = int_tick_label
    )
    
    plot!(p, x1_scaled, y1, label = legend1, lw = 2, marker = :circle)
    plot!(p, x2_scaled, y2, label = legend2, lw = 2, marker = :square)

    use_data_ticks = isnothing(ticks_from_data) ? (xscale != 1) : ticks_from_data
    if use_data_ticks
        ticks = sort(unique(vcat(x1_scaled, x2_scaled)))
        labels = int_tick_label.(ticks)
        plot!(p; xticks = (ticks, labels), xlims = (minimum(ticks), maximum(ticks)))
    end

    if !isnothing(tick_font)
        plot!(p; tickfont = font(tick_font))
    end
    
    # Determine output file name
    pdf_file = if isnothing(output_file)
        # Default: use file1 base name with _comparison.pdf extension
        base = splitext(file1)[1]
        "$(base)_comparison.pdf"
    else
        # Ensure .pdf extension
        endswith(output_file, ".pdf") ? output_file : "$(output_file).pdf"
    end
    
    # Ensure directory exists and save
    dir = dirname(pdf_file)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    savefig(p, pdf_file)
    println("Saved comparison plot to: $pdf_file")
    return p
end

"""
    plot_comparison_from_file(filename, legend1, legend2, xlabel, ylabel;
                              output_file=nothing, xscale=1.0,
                              ticks_from_data=nothing, tick_font=nothing,
                              legend_position=:best,
                              rotate_x_labels::Bool=false)

Reads columns 1 (x), 2 (y₁), and 3 (y₂) from a single file and plots the two
series together.

- `xscale`: multiply the x column by this factor before plotting.
- `ticks_from_data`: when `true`, place x ticks at every scaled data point.
- `legend_position`: legend placement symbol (default `:best`).
- `rotate_x_labels`: when `true`, rotate x-axis tick labels by 90°.
"""
function plot_comparison_from_file(filename::String, legend1::String, legend2::String,
                                   x_label::String, y_label::String;
                                   output_file=nothing, xscale::Real=2.0,
                                   ticks_from_data::Union{Bool,Nothing}=nothing,
                                   tick_font::Union{Nothing,Real}=nothing,
                                   legend_position::Symbol = :best,
                                   rotate_x_labels::Bool=false)
    x, y1, y2 = read_xyy_columns(filename)
    x_scaled = xscale == 1 ? x : x .* xscale

    p = plot(
        xlabel = x_label,
        ylabel = y_label,
        legend = legend_position,
        xrotation = rotate_x_labels ? 90 : 0,
        xformatter = int_tick_label
    )
    
    plot!(p, x_scaled, y1, label = legend1, lw = 2, marker = :circle)
    plot!(p, x_scaled, y2, label = legend2, lw = 2, marker = :square)

    use_data_ticks = isnothing(ticks_from_data) ? (xscale != 1) : ticks_from_data
    if use_data_ticks
        ticks = sort(unique(x_scaled))
        labels = int_tick_label.(ticks)
        plot!(p; xticks = (ticks, labels), xlims = (minimum(ticks), maximum(ticks)))
    end

    if !isnothing(tick_font)
        plot!(p; tickfont = font(tick_font))
    end

    # Determine output file name
    pdf_file = if isnothing(output_file)
        # Default: use input filename base with _comparison.pdf extension
        base = splitext(filename)[1]
        "$(base)_comparison.pdf"
    else
        # Ensure .pdf extension
        endswith(output_file, ".pdf") ? output_file : "$(output_file).pdf"
    end
    
    # Ensure directory exists and save
    dir = dirname(pdf_file)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
    savefig(p, pdf_file)
    println("Saved comparison plot to: $pdf_file")
    return p
end