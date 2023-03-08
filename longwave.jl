# Black body radiation

function planck(T, params)
    params.stephan*(T^4)
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
