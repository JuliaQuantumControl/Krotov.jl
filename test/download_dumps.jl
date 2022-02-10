using DrWatson
@quickactivate "KrotovTests"
import Downloads

DOWNLOADS = Dict(
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/DissGateOCT%23J_T%3DJ_T_re%23iter_stop%3D3000%23method%3Dkrotov.jld2" =>
        joinpath(datadir(), "DissGateOCT#J_T=J_T_re#iter_stop=3000#method=krotov.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_logisticsq.jld2" =>
        joinpath(datadir(), "parametrization#opt_result_logisticsq.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_positive.jld2" =>
        joinpath(datadir(), "parametrization#opt_result_positive.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_tanh.jld2" =>
        joinpath(datadir(), "parametrization#opt_result_tanh.jld2"),
    "https://github.com/JuliaQuantumControl/Krotov.jl/raw/data-dump/parametrization%23opt_result_tanhsq.jld2" =>
        joinpath(datadir(), "parametrization#opt_result_tanhsq.jld2"),
)

function download_dump(url, destination; force=false, verbose=true)
    if !isfile(destination) || force
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
