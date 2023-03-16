
function radiative_balance(params, t, temp)
    flux_down_sw = sw_down(params, t)
    flux_up_sw = sw_up(params, flux_down_sw)

    flux_tot_sw = flux_up_sw - flux_down_sw

    flux_down_lw = lw_down(temp, params)
    t_s = ((flux_down_lw[1]-flux_tot_sw[1]/params.emissiv)/params.stephan)^(1/4)
    flux_up_lw = lw_up(flux_down_lw, t_s, temp, params)

    flux_tot_lw = flux_up_lw - flux_down_lw

    flux_tot_uv = -uv_down(params, t)

    dT_sw = zeros(params.n)
    dT_lw = zeros(params.n)
    dT_uv = zeros(params.n)

    for i in 1:params.n
        dT_sw[i] = params.g*(flux_tot_sw[i] - flux_tot_sw[i+1])/(params.Cp*(params.p_int[i]-params.p_int[i+1]))
        dT_lw[i] = params.g*(flux_tot_lw[i] - flux_tot_lw[i+1])/(params.Cp*(params.p_int[i]-params.p_int[i+1]))
        dT_uv[i] = params.g*(flux_tot_uv[i] - flux_tot_uv[i+1])/(params.Cp*(params.p_int[i]-params.p_int[i+1]))
    end

    dT = dT_sw + dT_lw + dT_uv

    for i in 1:params.n
        temp[i] += params.dt*dT[i]
    end

    return t_s
end

function init_radiative(params)
    c_sw = -.5*log(params.coefvis)/(params.ps)
    c_lw = -log(params.coefir)/sqrt((params.ps)^2/(2*params.g))
    c_uv = -log(params.coefuv)/sqrt((params.ps)^2/(2*params.g))
    return (; params..., c_lw, c_sw, c_uv)
end
