using Distributions

struct ChangeFinder{T<:Real}
    discounting::T
    order::T
    smooth::Integer

    function ChangeFinder(discounting, order, smooth)
        if 0.0 < discounting < 1.0
            error("invalid discounting value. it must be 0.0 < discounting < 1.0.")
        end
        if order <= 0
            error("invalid order value. it must be order > 0.")
        end
        if smooth <= 1
            error("invalid smooth value. it must be smooth > 1.")
        end
        new(discounting, order, smooth)
    end
end

ChangeFinder(discounting::T, order::T, smooth::Integer) = ChangeFinder{T}(discounting, order, smooth)

function gaussian_likelihood(μ::Float64, σ::Float64, xt::Vector{Float64}, w::Vector{Float64})
    #m = MvNormal(μ::, σ::Float64)
    #pdf()
end

end
