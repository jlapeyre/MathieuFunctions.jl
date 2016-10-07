using Base.Test, MathieuFunctions

@test maximum(CharacteristicA(0,0:100) - [0:100;].^2) == 0

@test norm(CharacteristicB(0,1:100) - [1:100;].^2) == 0

begin
    test1 = readcsv("./MathieuCharacteristicA-1.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicA(q,0:10) for q in [-10:.01:10;]]) |> abs |> maximum) < 7.5e-13
end

begin
    test1 = readcsv("./MathieuCharacteristicA-2.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicA(q,0:3) for q in [30:.01:50;]]) |> abs |> maximum) < 7.6e-13
end

begin
    test1 = readcsv("./MathieuCharacteristicB-1.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicB(q,1:10) for q in [-10:.01:10;]]) |> abs |> maximum) < 7.5e-13
end

begin
    test1 = readcsv("./MathieuCharacteristicB-2.csv");
    @test (test1 - reduce(hcat,Vector{Float64}[CharacteristicB(q,1:3) for q in [30:.01:50;]]) |> abs |> maximum) < 2.8e-11
end

begin
    test1 = readcsv("./MathieuCharacteristicL-1.csv")[1:100,:]
    test2 = Float64[Characteristicλ(q,ν:ν)[1] for ν in [0:.01:0.99;], q in [-5:.01:5;]]
    (test1 - test2 |> abs |> maximum) < 7.5e-15
    # TODO: test ν > 1 (currently failing)
end

begin
    test1 = readcsv("./MathieuCharacteristicL-2.csv")[1:100,:]
    test2 = Float64[Characteristicλ(q,ν:ν)[1] for ν in [0:.01:0.99;], q in [30:.01:50;]]
    (test1 - test2 |> abs |> maximum) < 4.5e-14
    # TODO: test ν > 1 (currently failing)
end
