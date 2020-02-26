using SchistoIndividual, Test

@testset "make_death_rate_array" begin
    @test make_death_rate_array([6.56, 0.93, 0.3, 0.23, 0.27, 0.38, 0.44, 0.48,0.53, 0.65, 
    0.88, 1.06, 1.44, 2.1, 3.33, 5.29, 8.51, 13.66, 
    21.83, 29.98, 36.98], 1)[1] == (6.56/(1000*365))
end




@testset "find_death_rate" begin
    @test find_death_rate(0, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]) == 1
end


@testset "find_death_rate2" begin
    @test find_death_rate(1000, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]) == 18
    
end


@testset "find_death_rate3" begin
    @test find_death_rate(6, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]) == 3
    #x = "blah"
    #@test_throws MethodError find_death_rate(0, [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18])
end

@testset "make_age_contact_rate_array(max_age)" begin
    @test make_age_contact_rate_array(100)[1] == 0.22
end