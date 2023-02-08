using(Plots)


# Constants

params = (solarc = 1340, alb = 0.32, ps=101325, coefvis=0.99, coefir = 0.1, eps = 0.9, cpp = 1.0005, g = 9.81, stephan = 5.67e-8)
n = 10
c_sw = -.5*log(params.coefvis)/params.ps
c_lw = -log(params.coefir)/sqrt((params.ps)^2/(2*params.g))


# Pressure

p_int = zeros(n+1)

for i in 1:n+1
    p_int[i]=params.ps*(1-(i/(n+1)))
end

# Black body radiation

function planck(T)
    params.stephan*(T)^4
end

# Shortwave

flux_down_sw = zeros(n+1)
flux_up_sw = zeros(n+1)
tau_0_sw = zeros(n+1)
tau_inf_sw = zeros(n+1)

for i in 1:n+1
    tau_inf_sw[i] = exp(-c_sw*p_int[i])
    flux_down_sw[i] = tau_inf_sw[i]*params.solarc
    tau_0_sw[i] = exp(-(c_sw*params.ps+c_sw*p_int[i]))
    flux_up_sw[i] = params.alb*tau_0_sw[i]*flux_down_sw[1]
end

flux_tot_sw = flux_up_sw - flux_down_sw
    
# Longwave

function tau(i_1, i_2)
    exp(-sqrt(abs((c_lw)^2*((p_int[i_2])^2-(p_int[i_1])^2))))
end

function lw_down(flux_down_lw, temp)
    for i in 1:n+1
        for k in n:-1:i
            flux_down_lw[i] += (tau(k,i) - tau(k+1, i))*planck(temp[k])
        end
    end
    return flux_down_lw
end

function lw_up(flux_down , t_s, temp, flux_up_lw=zeros(n+1))
    flux_up_lw[1] = params.eps*params.stephan*t_s^4 + (1-params.eps)*flux_down[1]
    for i in 2:n+1
        for k in 1:i-1
            flux_up_lw[i] += (tau(i, k+1) - tau(i,k))*planck(temp[k]) + tau(1,i)*flux_up_lw[1]
        end
    end
    return flux_up_lw
end


# Changes in temperature 

function temp_ev(t_f,t=0, temp = 293*ones(n))
    dt = 10
    #Y = zeros(Int(ceil((t_f/dt)))+1, n+1)
    #Y[1,:] = vcat([293], temp)
    i = 2
    while t < t_f

        flux_down_lw = lw_down(zeros(n+1), temp)
        t_s = sqrt(sqrt(flux_down_lw[1]/params.stephan))
        flux_up_lw = lw_up(flux_down_lw, t_s, temp)

        flux_tot_lw = flux_up_lw - flux_down_lw

        dT_sw = zeros(n)
        dT_lw = zeros(n)

        for i in 1:n
            dT_sw[i] = params.g*(flux_tot_sw[i] - flux_tot_sw[i+1])/(params.cpp*(p_int[i]-p_int[i+1]))
            dT_lw[i] = params.g*(flux_tot_lw[i+1]-flux_tot_lw[i])/(params.cpp*(p_int[i+1]-p_int[i]))
        end
        dT = dT_sw + dT_lw
        for i in 1:n
            temp[i] += dt*dT[i]
        end
        #Y[i,:] = vcat([t_s], temp)
        #i += 1
        t += dt
        
    end
    println(temp)
    
    # plot(range(1,n+1), [Y[k,:] for k in 1:n+1])
end