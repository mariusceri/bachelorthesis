using(Plots)

# Black body radiation

function planck(T, params)
    params.stephan*(T^4)
end

# Shortwave

function sw_down(params)
    flux_down_sw = zeros(params.n +1)
    tau_inf_sw = zeros(params.n +1)
    for i in 1:params.n +1
        tau_inf_sw[i] = exp(-params.c_sw*params.p_int[i]/params.mu)
        flux_down_sw[i] = tau_inf_sw[i]*params.solarc*params.mu
    end
    return flux_down_sw
end

function sw_up(params, flux_down)
    flux_up_sw = zeros(params.n +1)
    tau_0_sw = zeros(params.n +1)
    flux_up_sw[1] = params.alb*flux_down[1]
    for i in 2:params.n+1
        tau_0_sw[i] = exp(-params.c_sw*(params.ps - params.p_int[i])/(3/5))
        flux_up_sw[i] = tau_0_sw[i]*flux_up_sw[1]
    end
    return flux_up_sw
end


# Longwave

function tau_lw(i_1, i_2, params)
    zdup = ((params.p_int[i_2])^2-(params.p_int[i_1])^2)/(2*params.g)
    return exp(-params.c_lw*sqrt(zdup))
end

function lw_down(temp, params)
    flux_down_lw = zeros(params.n +1)
    for i in 1:params.n +1  
        for k in params.n:-1:i   # k>=i 
            flux_down_lw[i] += (tau_lw(k,i, params) - tau_lw(k+1, i, params))*planck(temp[k], params)
        end
    end
    return flux_down_lw
end

function lw_up(flux_down , t_s, temp, params)
    flux_up_lw = zeros(params.n +1)
    flux_up_lw[1] = params.emissiv*planck(t_s, params) + (1-params.emissiv)*flux_down[1]
    for i in 2:params.n +1
        flux_up_lw[i] = tau_lw(i, 1, params)*flux_up_lw[1]

        for k in 1:i-1 # i>=k+1
            flux_up_lw[i] += (tau_lw(i, k+1, params) - tau_lw(i,k, params))*planck(temp[k], params)
        end
    end
    return flux_up_lw
end

# Changes in temperature 

function temp_ev(t_f, params)
    temp = 293*ones(params.n)
    dt = 3600
    c = 10                               # number of curves to plot
    Y = zeros(c, params.n +1)            # initialisation of the vector for the plot   
    i = 1                                # iterator for the plot 
    l = [elem for elem in 1:t_f/c:t_f]

    for t in 1:dt:t_f

        flux_down_sw = sw_down(params)
        flux_up_sw = sw_up(params, flux_down_sw)

        flux_tot_sw = flux_up_sw - flux_down_sw

        flux_down_lw = lw_down(temp, params)
        t_s = ((flux_down_lw[1]-flux_tot_sw[1]/params.emissiv)/params.stephan)^(1/4)
        flux_up_lw = lw_up(flux_down_lw, t_s, temp, params)

        flux_tot_lw = flux_up_lw - flux_down_lw

        dT_sw = zeros(params.n)
        dT_lw = zeros(params.n)

        for i in 1:params.n
            dT_sw[i] = params.g*(flux_tot_sw[i] - flux_tot_sw[i+1])/(params.cpp*(params.p_int[i]-params.p_int[i+1]))
            dT_lw[i] = params.g*(flux_tot_lw[i] - flux_tot_lw[i+1])/(params.cpp*(params.p_int[i]-params.p_int[i+1]))
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
    params = (n = 50, solarc = 1340/4, alb = 0.32, ps=101325, coefvis=0.9, coefir = 0.08, mu = 0.7,  emissiv = 0.9, cpp = 1005, g = 9.81, stephan = 5.67e-8)
    c_sw = -.5*log(params.coefvis)/(params.ps)
    c_lw = -log(params.coefir)/sqrt((params.ps)^2/(2*params.g))
    params = (; params..., c_lw, c_sw)

    # Pressure
    p_int = zeros(params.n+1)

    for i in 1:params.n +1
        p_int[i]=params.ps*(1-(i/(params.n+1)))
    end

    params = (; params..., p_int)

    temp_ev(3600*15000, params)
end

main()