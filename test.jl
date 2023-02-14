include("temp.jl")

# Shortwave
flux_down_sw = zeros(params.n +1)
flux_up_sw = zeros(params.n +1)
tau_0_sw = zeros(params.n +1)
tau_inf_sw = zeros(params.n +1)

for i in 1:params.n +1
    tau_inf_sw[i] = exp(-c_sw*p_int[i])
    flux_down_sw[i] = tau_inf_sw[i]*params.solarc
    tau_0_sw[i] = exp(-(c_sw*params.ps+c_sw*p_int[i]))
    flux_up_sw[i] = params.alb*tau_0_sw[i]*flux_down_sw[1]
end


println("Flux down shortwave:")
println(flux_down_sw)
println("Flux up shortwave:")
println(flux_up_sw)

println("Flux down longwave:")
test_flux_down = lw_down(zeros(params.n +1), 300*ones(params.n))
println("Flux up longwave:")
test_flux_up = lw_up(flux_down , 300, 300*ones(params.n))