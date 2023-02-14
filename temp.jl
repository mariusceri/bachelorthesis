using(Plots)

# Black body radiation

function planck(T)
    params.stephan*(T)^4
end

# Longwave

function tau_lw(i_1, i_2)
    exp(-sqrt(abs((params.c_lw)^2*((p_int[i_2])^2-(p_int[i_1])^2))))
end

function lw_down(flux_down_lw, temp)
    for i in 1:params.n +1
        for k in params.n:-1:i
            flux_down_lw[i] += (tau_lw(k,i) - tau_lw(k+1, i))*planck(temp[k])
        end
    end
    println(flux_down_lw)
    return flux_down_lw
end

function lw_up(flux_down , t_s, temp, flux_up_lw=zeros(params.n +1))
    flux_up_lw[1] = params.eps*params.stephan*t_s^4 + (1-params.eps)*flux_down[1]
    for i in 2:params.n +1
        for k in 1:i-1
            flux_up_lw[i] += (tau_lw(i, k+1) - tau_lw(i,k))*planck(temp[k]) + tau_lw(1,i)*flux_up_lw[1]
        end
    end
    println(flux_up_lw)
    return flux_up_lw
end

# Changes in temperature 

function temp_ev(t_f, temp = 293*ones(params.n))
    dt = 40
    c = 10                               # number of curves to plot
    Y = zeros(c, params.n +1)            # initialisation of the vector for the plot   
    i = 1                                # iterator for the plot 
    l = [elem for elem in 1:t_f/c:t_f]

    for t in 1:dt:t_f

        flux_down_lw = lw_down(zeros(params.n +1), temp)
        t_s = sqrt(sqrt(flux_down_lw[1]/params.stephan))
        flux_up_lw = lw_up(flux_down_lw, t_s, temp)

        flux_tot_lw = flux_up_lw - flux_down_lw

        dT_sw = zeros(params.n)
        dT_lw = zeros(params.n)

        for i in 1:params.n
            dT_sw[i] = params.g*(flux_tot_sw[i+1] - flux_tot_sw[i])/(params.cpp*(p_int[i+1]-p_int[i]))
            dT_lw[i] = params.g*(flux_tot_lw[i+1]-flux_tot_lw[i])/(params.cpp*(p_int[i+1]-p_int[i]))
        end

        dT = dT_sw + dT_lw

        for i in 1:params.n
            temp[i] += dt*dT[i]
        end

        if t in l
            Y[i,:] = vcat([t_s], temp)
            i += 1
        end

    end
    println(temp)
    plot([Y[k,:] for k in 1:c], range(1,params.n +1))
end

# Execution

function main()
    params = (n = 5, solarc = 1340, alb = 0.32, ps=101325, coefvis=0.9, coefir = 0.01, eps = 0.9, cpp = 1.0005, g = 9.81, stephan = 5.67e-8)
    c_sw = -.5*log(params.coefvis)/params.ps
    c_lw = -log(params.coefir)/sqrt((params.ps)^2/(2*params.g))
    params = (; params..., c_lw, c_sw)

    # Pressure
    p_int = zeros(params.n+1)

    for i in 1:params.n +1
        p_int[i]=params.ps*(1-(i/(params.n+1)))
    end

    println(p_int)

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
end