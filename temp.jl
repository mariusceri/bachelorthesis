using(Plots)
using(LaTeXStrings)

include("shortwave.jl")
include("longwave.jl")
include("radiative_balance.jl")
include("thermodynamics.jl")
include("dry_adjustment.jl")

# Changes in temperature 

function temp_ev(t_f, params, temp ;
                c=10)                    # number of curves to plot
    Y = zeros(c+1, params.n +1)            # initialisation of the vector for the plot   
    i = 1                                # iterator for the plot 
    l = [elem for elem in 1:t_f/c:t_f]

    for t in 1:params.dt:t_f

        t_s = radiative_balance(params, t, temp)
        params.adjust && adjust_Nlayers!(params, params.p_layer, temp)

        if t in l
            Y[i,:] = vcat([t_s], temp)
            i += 1
        end
    Y[i,:] = vcat([t_s], temp)
    end
    return temp, Y
end

# Execution

positive(x) = (x+abs(x))/2
diurnal_cycle(t) = positive(cos(t/24/3600*2pi))
no_diurnal_cycle(t)=1

function main()
    params = (n = 100, dt=3600, # numerical parameters
            solarc = 1340/4, diurnal_cycle = no_diurnal_cycle,
            ps = 101325, g = 9.81,
            Pr = 101325.0, Cp = 1005, R = 287, # perfect gas
            coefvis=0.9, coefuv=1e-6, alb = 0.32, mu = 0.7,  # parameters for SW
            stephan = 5.67e-8, coefir = 0.2, emissiv = 0.9, # parameters for LW
            adjust=false,
            )
    params = init_radiative(params)

    # Pressure
    p_int   = zeros(params.n+1)
    p_layer = zeros(params.n)

    for i in 1:params.n +1
        p_int[i]=params.ps*(1-(i/(params.n+1)))
    end
    for i in 1:params.n
        p_layer[i] = (p_int[i]+p_int[i+1])/2
    end

    params = (; params..., p_int, p_layer)

    temp = 293*ones(params.n)
    temp, Y = temp_ev(3600*24*500, params, temp)
    #temp, Y = temp_ev(3600*24*2, params, temp, c=24*2)

    println(temp)
    display(plot([Y[k,:] for k in 1:11], range(1,params.n +1), label =[L"t=0" L"t=t_f/10" L"t=2t_f/10" L"t=3t_f/10" L"t=4t_f/10" L"t=5t_f/10" L"t=6t_f/10" L"t=7t_f/10" L"t=8t_f/10" L"t=9t_f/10" L"t=t_f"]))
    #        label = [L"T(0)" L"T(\frac{t_{f}}{10})" L"T(\frac{2t_{f}}{10})" 
    #        L"T(\frac{3t_{f}}{10})" L"T(\frac{4t_{f}}{10})" L"T(\frac{5t_{f}}{10})" 
    #        L"T(\frac{6t_{f}}{10})" L"T(\frac{7t_{f}}{10})" L"T(\frac{8t_{f}}{10})" L"T(\frac{9t_{f}}{10})" L"T(t_f)"]))
    #nn=60
    #display(plot(Y[:,1]))

    #display(plot([potential_temperature(params, params.p_layer[1:nn], Y[k,2:nn+1]) for k in axes(Y,1)], params.p_layer[1:nn] ; yflip=true))
    xlabel!("Temperature (K)")
    ylabel!("Altitude (number of layers)")
end

main()

#profview