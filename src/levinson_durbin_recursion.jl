module LevinsonDurbin

function solve(toeplitz_eles::Vector{Float64}) :: Vector{Float64}
    solutions, extra_ele = solve_size_one(toeplitz_eles)
    lpc_dim = size(toeplitz_eles, 1) - 2
    final_solutions, a = solve_recursively(toeplitz_eles, solutions, extra_ele, lpc_dim)
    return final_solutions
end

function solve_size_one(toeplitz_eles::Vector{Float64}) :: Tuple{Vector{Float64}, Float64}
    solutions = [1.0, - toeplitz_eles[2] / toeplitz_eles[1]]
    extra_element = toeplitz_eles[1] + toeplitz_eles[2] * solutions[2]
    return solutions, extra_element
end

function solve_recursively(
    toeplitz_eles::Vector{Float64},
    initial_solutions::Vector{Float64},
    initial_extra_ele::Float64,
     lpc_dim::Int64) :: Tuple{Vector{Float64}, Float64}

    extra_element = initial_extra_ele
    solutions = initial_solutions
    for k = 2:lpc_dim
        lambda = calc_lambda(k, toeplitz_eles, solutions, extra_element)
        extended_solution = push!(solutions, 0.0) # extend the solution.
        r_extended_solution = reverse(extended_solution)

        solutions = extended_solution + lambda * r_extended_solution
        extra_ele = (1.0 - lambda^2) * extra_ele # next extra element.
    end
    return solutions, extra_element
end

function calc_lambda(k::Int64, toeplitz_eles::Vector{Float64}, solutions::Vector{Float64}, extra_ele::Float64) :: Float64
    tmp_value = 0.0
    for j = 1:k
        tmp_value += (- solutions[j] * toeplitz_eles[k + 2 - j])
    end
    return tmp_value / extra_ele
end

end
