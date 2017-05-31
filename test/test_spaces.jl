# test_spaces.jl

function test_spaces()
    test_basic_spaces()
    test_euclidean_spaces()
end

function test_basic_spaces()
end

function test_euclidean_spaces()
end

function test_isomorphism(space1, space2)
    @test isomorphic(space1, space2)
end

function test_embedding(space1, space2)
    @test embedded(space1, space2)
end
