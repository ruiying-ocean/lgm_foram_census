#!/usr/bin/env julia
#
# DIVA interpolation of planktonic foraminifera census data onto a 1×1° global grid.
# Produces CF-compliant NetCDF files for LGM and PI (ForCenS) time slices.
#
# Usage:  julia code/diva_interpolation.jl

using CSV, DataFrames, Statistics, Downloads, Dates
using DIVAnd, NCDatasets

# ── Configuration ────────────────────────────────────────────────────────────
const GRID_LON = -180:1:179        # 360 cells
const GRID_LAT = -90:1:90          # 181 cells
const CORR_LEN = 1000e3            # correlation length [m]
const EPSILON2 = 0.1               # noise-to-signal (ε² = 1/SNR)
const MIN_OBS  = 10                # skip species with fewer valid observations

const BATHY_URL  = "https://dox.ulg.ac.be/index.php/s/U0pqyXhcQrXjEUX/download"
const BATHY_FILE = joinpath(@__DIR__, "..", "data", "gebco_30sec_16.nc")

const LGM_CSV     = joinpath(@__DIR__, "..", "tidy", "lgm_sp_r_tidy.csv")
const FORCENS_CSV  = joinpath(@__DIR__, "..", "tidy", "forcens_sp_r_tidy.csv")
const LGM_NC      = joinpath(@__DIR__, "..", "tidy", "lgm_sp_r_diva.nc")
const FORCENS_NC   = joinpath(@__DIR__, "..", "tidy", "forcens_sp_r_diva.nc")

# ── Bathymetry & mask ────────────────────────────────────────────────────────

function download_bathymetry()
    if !isfile(BATHY_FILE)
        @info "Downloading GEBCO bathymetry → $BATHY_FILE"
        Downloads.download(BATHY_URL, BATHY_FILE)
    end
end

function setup_grid()
    xi, yi = ndgrid(collect(Float64, GRID_LON), collect(Float64, GRID_LAT))
    pm, pn = DIVAnd_metric(xi, yi)
    return xi, yi, pm, pn
end

function create_mask()
    _, _, mask = load_mask(BATHY_FILE, true, GRID_LON, GRID_LAT, 0)
    return mask
end

# ── Data loading ─────────────────────────────────────────────────────────────

function load_lgm_data()
    df = CSV.read(LGM_CSV, DataFrame; missingstring="")

    # Species columns: columns 6 onward (after 5 metadata columns)
    sp_names = names(df)[6:end]

    # Fix lon > 180
    df.Longitude = ifelse.(df.Longitude .> 180, df.Longitude .- 360, df.Longitude)

    # Average per core (Event): collapse depth replicates
    # Use safe mean that returns missing when all values are missing
    safemean(x) = let v = collect(skipmissing(x)); isempty(v) ? missing : mean(v) end
    gd = groupby(df, :Event)
    agg = combine(gd,
        :Latitude  => first => :Latitude,
        :Longitude => first => :Longitude,
        [sp => safemean => sp for sp in sp_names]...
    )

    lon = Float64.(agg.Longitude)
    lat = Float64.(agg.Latitude)
    sp_data = Dict{String,Vector{Union{Missing,Float64}}}()
    for sp in sp_names
        sp_data[sp] = Vector{Union{Missing,Float64}}(agg[!, sp])
    end

    @info "LGM: $(length(lon)) cores, $(length(sp_names)) species"
    return lon, lat, sp_names, sp_data
end

function load_forcens_data()
    df = CSV.read(FORCENS_CSV, DataFrame; missingstring="")

    # Species columns start at column 22
    sp_names = names(df)[22:end]

    lon = Float64.(df.Longitude)
    lat = Float64.(df.Latitude)

    # Fix any lon > 180
    lon .= ifelse.(lon .> 180, lon .- 360, lon)

    sp_data = Dict{String,Vector{Union{Missing,Float64}}}()
    for sp in sp_names
        sp_data[sp] = Vector{Union{Missing,Float64}}(df[!, sp])
    end

    @info "ForCenS: $(length(lon)) samples, $(length(sp_names)) species"
    return lon, lat, sp_names, sp_data
end

# ── Interpolation ────────────────────────────────────────────────────────────

function interpolate_species(xi, yi, pm, pn, mask,
                             lon, lat, values::Vector{Union{Missing,Float64}};
                             corr_len=CORR_LEN, epsilon2=EPSILON2)

    # Keep only non-missing observations
    valid = .!ismissing.(values)
    nvalid = count(valid)

    if nvalid < MIN_OBS
        return fill(NaN, size(xi))
    end

    obs_lon = lon[valid]
    obs_lat = lat[valid]
    obs_val = Float64.(values[valid])

    fi, _ = DIVAndrun(
        mask, (pm, pn), (xi, yi),
        (obs_lon, obs_lat), obs_val,
        (corr_len, corr_len), epsilon2;
        moddim = [360.0, 0.0]
    )

    return fi
end

function renormalize!(fields::Dict{String,Matrix{Float64}}, mask)
    nx, ny = size(first(values(fields)))
    sp_keys = collect(keys(fields))

    for i in 1:nx, j in 1:ny
        if !mask[i, j]
            for sp in sp_keys
                fields[sp][i, j] = NaN
            end
            continue
        end

        # Clamp negatives
        for sp in sp_keys
            v = fields[sp][i, j]
            if isnan(v) || v < 0
                fields[sp][i, j] = 0.0
            end
        end

        total = sum(fields[sp][i, j] for sp in sp_keys)

        if total > 0
            for sp in sp_keys
                fields[sp][i, j] /= total
            end
        end
    end
end

# ── NetCDF output ────────────────────────────────────────────────────────────

function save_netcdf(filename, sp_names, fields, mask;
                     title="DIVAnd interpolated foraminifera")
    lon_vals = collect(Float64, GRID_LON)
    lat_vals = collect(Float64, GRID_LAT)

    isfile(filename) && rm(filename)

    NCDataset(filename, "c"; attrib=Dict(
        "title"       => title,
        "institution" => "DIVAnd interpolation",
        "source"      => "LGM foram census project",
        "history"     => "Created $(Dates.today()) by diva_interpolation.jl",
        "Conventions" => "CF-1.8"
    )) do ds
        defDim(ds, "lon", length(lon_vals))
        defDim(ds, "lat", length(lat_vals))

        lon_var = defVar(ds, "lon", Float64, ("lon",);
            attrib=Dict("units" => "degrees_east",
                        "long_name" => "longitude",
                        "standard_name" => "longitude"))
        lon_var[:] = lon_vals

        lat_var = defVar(ds, "lat", Float64, ("lat",);
            attrib=Dict("units" => "degrees_north",
                        "long_name" => "latitude",
                        "standard_name" => "latitude"))
        lat_var[:] = lat_vals

        for sp in sp_names
            # sanitize variable name: replace dots/spaces with underscores
            varname = replace(sp, r"[. ]" => "_")
            data = fields[sp]

            # Convert NaN to fill value for masked cells
            v = defVar(ds, varname, Float32, ("lon", "lat");
                fillvalue=Float32(-9999.0),
                attrib=Dict("long_name" => sp,
                            "units" => "1",
                            "comment" => "relative abundance (0-1)"))

            out = Array{Union{Missing,Float32}}(undef, size(data))
            for i in eachindex(data)
                out[i] = isnan(data[i]) ? missing : Float32(data[i])
            end
            v[:, :] = out
        end
    end
    @info "Saved $filename"
end

# ── Main ─────────────────────────────────────────────────────────────────────

function main()
    @info "Starting DIVA interpolation"

    download_bathymetry()
    xi, yi, pm, pn = setup_grid()
    mask = create_mask()
    @info "Grid: $(size(xi)), ocean cells: $(count(mask))"

    # ── LGM ──────────────────────────────────────────────────────────────
    lon_lgm, lat_lgm, sp_lgm, data_lgm = load_lgm_data()

    fields_lgm = Dict{String,Matrix{Float64}}()
    for (i, sp) in enumerate(sp_lgm)
        @info "  LGM [$i/$(length(sp_lgm))] $sp"
        fields_lgm[sp] = interpolate_species(
            xi, yi, pm, pn, mask,
            lon_lgm, lat_lgm, data_lgm[sp])
    end

    renormalize!(fields_lgm, mask)
    save_netcdf(LGM_NC, sp_lgm, fields_lgm, mask;
                title="LGM planktonic foraminifera relative abundance (DIVAnd)")

    # ── ForCenS (PI) ─────────────────────────────────────────────────────
    lon_fc, lat_fc, sp_fc, data_fc = load_forcens_data()

    fields_fc = Dict{String,Matrix{Float64}}()
    for (i, sp) in enumerate(sp_fc)
        @info "  ForCenS [$i/$(length(sp_fc))] $sp"
        fields_fc[sp] = interpolate_species(
            xi, yi, pm, pn, mask,
            lon_fc, lat_fc, data_fc[sp])
    end

    renormalize!(fields_fc, mask)
    save_netcdf(FORCENS_NC, sp_fc, fields_fc, mask;
                title="PI (ForCenS) planktonic foraminifera relative abundance (DIVAnd)")

    @info "Done."
end

main()
