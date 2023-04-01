import Downloads

const PROJECTDIR = dirname(Base.active_project())
projectdir(names...) = joinpath(PROJECTDIR, names...)
datadir(names...) = projectdir("data", names...)


DOWNLOADS = Dict(
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/DissGateOCT%23J_T%3DJ_T_re%23iter_stop%3D3000%23method%3Dkrotov.jld2" =>
        datadir("DissGateOCT#J_T=J_T_re#iter_stop=3000#method=krotov.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_logisticsq.jld2" =>
        datadir("parametrization#opt_result_logisticsq.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_positive.jld2" =>
        datadir("parametrization#opt_result_positive.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_tanh.jld2" =>
        datadir("parametrization#opt_result_tanh.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_tanhsq.jld2" =>
        datadir("parametrization#opt_result_tanhsq.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/PE_OCT.jld2" =>
        datadir("PE_OCT.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/PE_OCT_direct.jld2" =>
        datadir("PE_OCT_direct.jld2"),
)

function download_dump(url, destination; force=false, verbose=true)
    if !isfile(destination) || force
        mkpath(dirname(destination))
        verbose && (@info "Downloading $url => $destination")
        Downloads.download(url, destination)
    else
        verbose && (@info "$destination OK")
    end
end

@info "Download Dumps"
for (url, destination) in DOWNLOADS
    download_dump(url, destination)
end
