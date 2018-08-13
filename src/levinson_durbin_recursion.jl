
function levinson_durbin_recursion(toeplitz_elements::Vector{Float64})
    return toeplitz_elements
end

function solve_size_one(toeplitz_elements::Vector{Float64})
    solutions = [1.0, - toeplitz_elements[2] / toeplitz_elements[1]]
    extra_element = toeplitz_elements[1] + toeplitz_elements[2] * solutions[2]
    return solutions, extra_element
end

function solve_recursively(
    initial_solutions::Vector{Float64},
    initial_extra_ele::Float64, lpc_dim:Int64)

    extra_element = initial_extra_ele
    solutions = initial_solutions
    for k = 1:(lpc_dim-1)
        lambda_value = calculate_lambda(k, solutions, extra_element)
        extended_solution = extend_solution(solutions)
        r_extended_solution = reverse(extended_solution)

        solutions = extended_solution + lambda_value * r_extended_solution
        extra_element = calculate_extra_element(extra_element, lambda_value)
    end
    return solutions, extra_element
end

function extend_solutions(previous_soloution::Vector{Float64})
    return push!(previous_soloution, 0.0)
end

function calculate_extra_element(previous_extea_ele::Float64, lambda_value::Flout64)
    return (1.0 - lambda_value^2) * previous_extra_ele
end

function calculate_lambda(k::Int64, toeplitz_elements::Vector{Float64}, solutions::Vector{Float64}, extra_ele::Float64)
    tmp_value = 0.0
    for j = 1:k
        tmp_value += (- solutions[j] * toeplitz_elements[k + 1 - j])
    end
    return tmp_value / extra_element
end
