using ChangeFinder.LevinsonDurbin
using Test

@testset "solve" begin
    toeplitz_elements = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    expected_solution = [1.0, -1.2, 4.44089e-17, 3.88578e-16, -0.2]
    @test expected_solution ≈ LevinsonDurbin.solve(toeplitz_elements)
end
