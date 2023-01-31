# Constants

solarc = 1340
alb = 0.32
N = 20
ps = 101325
coefvis = 0.99
coefir = 0.99
eps = 0.92
cpp = 1.0005
g = 9.81
c_sw = -.5*log(coefvis)/ps
c_lw = -log(coefir)/sqrt(ps*ps/(2*g))
stephan = 5.67e-8


# Initialisation of arrays

flux_down_sw = zeros(N+1)
flux_up_sw = zeros(N+1)
tau_0_sw = zeros(N+1)
tau_inf_sw = zeros(N+1)
p_int = zeros(N+1)

temp = 293*ones(N)

# Pressure

for i in 1:N+1
    p_int[i]=ps*(1-(i/N))
end

# Black body radiation

function planck(T)
    stephan*(T)^4
end

# Computation of shortwave fluxes

for i in 1:N+1
    tau_inf_sw[i] = exp(-c_sw*p_int[i])
    flux_down_sw[i] = tau_inf_sw[i]*solarc
end 

for i in 1:N+1
    tau_0_sw[i] = exp(-(c_sw*ps+c_sw*p_int[i]))
    flux_up_sw[i] = alb*tau_0_sw[i]*flux_down_sw[1]
end

flux_tot_sw = flux_up_sw - flux_down_sw

function tau(i_1, i_2)
    exp(-sqrt(abs(c_lw*c_lw*(p_int[i_2]*p_int[i_2]-p_int[i_1]*p_int[i_1]))))
end

# Changes in temperature 

function temp_ev(t_f, t=0, temp = 293*ones(N))
    dt = 50
    while t < t_f
        flux_down_lw = zeros(N+1)
        flux_up_lw = zeros(N+1)
        dT_sw = zeros(N)
        dT_lw = zeros(N)
        for i in 1:N+1
            for k in i+1:N
                flux_down_lw[i] += (tau(i,k) - tau(i, k+1))*planck(temp[k]) 
            end 
        end 
        
        for i in 1:N+1
            for k in 1:i-1
                flux_up_lw[i] += (tau(k+1,i) - tau(k,i))*planck(temp[k]) + tau(1,i)*(eps*planck(temp[1]) + (1-eps)*flux_down_lw[1])
            end
        end

        flux_tot_lw = flux_up_lw - flux_down_lw

        for i in 1:N
            dT_sw[i] = g*(flux_tot_sw[i] - flux_tot_sw[i+1])/(cpp*(p_int[i]-p_int[i+1]))
            dT_lw[i] = g*(flux_tot_lw[i+1]-flux_tot_lw[i])/(cpp*(p_int[i+1]-p_int[i]))
        end
        dT = dT_sw + dT_lw
        for i in 1:N
            temp[i] += dt*dT[i]
        end
        t += dt
    end
    println(temp)
end