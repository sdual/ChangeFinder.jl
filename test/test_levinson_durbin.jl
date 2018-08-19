using ChangeFinder.LevinsonDurbin
using Test

@testset "solve" begin
    toeplitz_eles = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    expected = [1.0, -1.2, 4.44089e-17, 3.88578e-16, -0.2]
    @test expected ≈ LevinsonDurbin.solve(toeplitz_eles)
end

@testset "solve_size_one" begin
    toeplitz_eles = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    expected_solutions = [1.0, -2.0]
    expected_extra_ele = -3.0
    actual = LevinsonDurbin.solve_size_one(toeplitz_eles)
    @test expected_solutions ≈ actual[1]
    @test expected_extra_ele ≈ actual[2]
end

@testset "solve_recursively" begin
    toeplitz_eles = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    initial_solutions = [1.0, -2.0]
    initial_extra_ele = -3.0
    lpc_dim = 4
    expected_solutions = [1.0, -1.2, 4.44089e-17, 3.88578e-16, -0.2]
    expected_extra_ele = -2.3999999999999995
    actual = LevinsonDurbin.solve_recursively(toeplitz_eles, initial_solutions, initial_extra_ele, lpc_dim)
    @test expected_solutions ≈ actual[1]
    @test expected_extra_ele ≈ actual[2]
end

@testset "calc_lambda" begin
    k = 2
    toeplitz_eles = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    solutions = [1.0, -2.0]
    extra_ele = -3.0
    expected = -0.3333333333333333
    @test expected ≈ LevinsonDurbin.calc_lambda(k, toeplitz_eles, solutions, extra_ele)
end
