function adjust_Nlayers!(params, p, T)
    theta = potential_temperature.(Ref(params), p, T)
    for n in (2:params.n-1)
        if theta[1] > theta[n]
            coef = exner(params, p)
            theta[1:n] .= sum(T[1:n])/sum(coef[1:n])
            T[1:n] = inverse_PT(params, p, theta[1])[1:n]
        end
    end
end
