# Constants

solarc = 1340
alb = 0.32
N = 5
ps = 101325
coefir = 0.9
eps = 0.99
cpp = 1.0005
g = 9.81
c_lw = -log(coefir)/sqrt(ps*ps/(2*g))
stephan = 5.67e-8

# Initialisation of arrays

flux_down_lw = zeros(N+1)
flux_up_lw = zeros(N+1)
p_int = zeros(N+1)
dT_lw = zeros(N)
temp = 293*ones(N)
zplanck = zeros(N)

# Pressure

for i in 1:N+1
    p_int[i]=ps*(1-(i/N))
end

# Black body radiation

for i in 1:N
    zplanck[i] = stephan*(temp[i])^4
end

# Computation of fluxes 

function tau(i_1, i_2)
    exp(-sqrt(abs(c_lw*c_lw*(p_int[i_2]*p_int[i_2]-p_int[i_1]p_int[i_1]))))
end

for i in 1:N+1
    for k in i+1:N
        flux_down_lw[i] += (tau(i,k) - tau(i, k+1))*zplanck(k) 
    end 
end 

for i in 1:N+1
    for k in 1:i-1
        flux_up_lw[i] += (tau(k+1,i) - tau(k,i))*zplanck(k) + tau(1,i)(eps*zplanck[1] + (1-eps)*flux_down_lw[1])
    end

flux_tot_lw = flux_up_lw - flux_down_lw

# Change in temperature 

for i in 1:N
    dT_lw[i] = (g/cpp)*(flux_tot_lw[i] - flux_tot_lw[i+1])/(p_int[i]-p_int[i+1])
end

function temp_ev_lw(t_f, t=0, temp = 293*ones(N))
    dt = 1
    while t < t_f
        for i in 1:N
            temp[i] += dt*dT_lw[i]
            t += dt
        end
        t += dt
    end
    println(temp)
end