# Constants

solarc = 1340
alb = 0.32
N = 10
ps = 101325
coefvis = 0.9
cpp = 1.0005
g = 9.81
c_sw = -.5*log(coefvis)/ps

# Initialisation of arrays

flux_down_sw = zeros(N+1)
flux_up_sw = zeros(N+1)
tau_0_sw = zeros(N+1)
tau_inf_sw = zeros(N+1)
p_int = zeros(N+1)
dT_sw = zeros(N)


# Pressure

for i in 1:N+1
    p_int[i]=ps*(1-(i/N))
end

# Computation of fluxes 

for i in 1:N+1
    tau_inf_sw[i] = exp(-c_sw*p_int[i])
    flux_down_sw[i] = tau_inf_sw[i]*solarc
end 

for i in 1:N+1
    tau_0_sw[i] = exp(-(c_sw*ps+c_sw*p_int[i]))
    flux_up_sw[i] = alb*tau_0_sw[i]*flux_down_sw[1]
end

flux_tot_sw = flux_up_sw - flux_down_sw

# Change in temperature 

for i in 1:N
    dT_sw[i] = (g/cpp)*(flux_tot_sw[i] - flux_tot_sw[i+1])/(p_int[i]-p_int[i+1])
end

function temp_ev_sw(t_f, t=0, temp_sw = 293*ones(N))
    dt = 1
    while t < t_f
        for i in 1:N
            temp_sw[i] += dt*dT_sw[i]
            t += dt
        end
        t += dt
    end
    println(temp_sw)
end