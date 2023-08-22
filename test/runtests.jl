using Test, MathieuFunctions,MathieuF
using LinearAlgebra
using DelimitedFiles

readcsv(f) = DelimitedFiles.readdlm(f, ',')

function tapprox(a, b; atol=1e-15)
    normval = norm(a - b, Inf)
    @info "normval = $normval"
    isapprox(a, b; norm= x -> norm(x, Inf), atol = atol)
end

@testset "basic" begin
    @test maximum([MathieuCharA(ν,0) for ν in 0:100] - [0:100;].^2) == 0
    @test norm([MathieuCharB(ν,0) for ν in 1:100] - [1:100;].^2) == 0
end


filename = "MathieuCharacteristicA-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r=[MathieuCharA(ν,q) for ν in 0:10, q in -10:.01:10]
    @test tapprox(test1, r;  atol=7.5e-13)
end

filename = "MathieuCharacteristicA-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = [MathieuCharA(ν,q) for ν in 0:3, q in 30:.01:50]
    @test tapprox(test1, r;  atol=7.6e-13) # NOTE: was 7.5e-13
end

filename = "MathieuCharacteristicB-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = [MathieuCharB(ν,q) for ν in 1:10, q in -10:.01:10]
    @test tapprox(test1, r;  atol=7.5e-13)
end

filename = "MathieuCharacteristicB-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    r = [MathieuCharB(ν,q) for ν in 1:3, q in 30:.01:50]
    @test tapprox(test1, r; atol=2.8e-11)
end

filename = "MathieuCharacteristicL-1.csv"
@testset "$filename" begin
    test1 = readcsv(filename)[1:100,:]
    test2 = Float64[MathieuCharA(ν,q) for ν in [0:.01:0.99;], q in -5:.01:5]
    @test tapprox(test1, test2, atol=7.5e-15)
end

filename = "MathieuCharacteristicL-2.csv"
@testset "$filename" begin
    test1 = readcsv(filename)[1:100,:]
    test2 = Float64[MathieuCharA(ν,q) for ν in [0:.01:0.99;], q in 30:.01:50]
    @test tapprox(test1, test2, atol=4.5e-14)
end

filename = "MathieuCharacteristicL-3.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    test2 = Float64[MathieuCharA(ν,q) for ν in [1.01:.01:1.99;], q in -5:.01:5]
    test3 = Float64[MathieuCharB(ν,q) for ν in [1.01:.01:1.99;], q in -5:.01:5]
    @test tapprox(test1, test2, atol=6e-14)
    @test tapprox(test1, test3, atol=6e-14)
end

filename = "MathieuCharacteristicL-4.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    test2 = Float64[MathieuCharA(ν,q) for ν in [1.01:.01:1.99;], q in 30:.01:50]
    test3 = Float64[MathieuCharB(ν,q) for ν in [1.01:.01:1.99;], q in 30:.01:50]
    @test tapprox(test1, test2, atol=6e-14)
    @test tapprox(test1, test3, atol=6e-14)
end

filename = "MathieuCharacteristicL-5.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    test2 = Float64[MathieuCharA(ν,q) for ν in [20.01:.01:20.99;], q in -5:.01:5]
    test3 = Float64[MathieuCharB(ν,q) for ν in [20.01:.01:20.99;], q in -5:.01:5]
    @test tapprox(test1, test2, atol=6e-13)
    @test tapprox(test1, test3, atol=6e-13)
end

filename = "MathieuCharacteristicL-6.csv"
@testset "$filename" begin
    test1 = readcsv(filename)
    test2 = Float64[MathieuCharA(ν,q) for ν in [20.01:.01:20.99;], q in 30:.01:50]
    test3 = Float64[MathieuCharB(ν,q) for ν in [20.01:.01:20.99;], q in 30:.01:50]
    @test tapprox(test1, test2, atol=7e-13)
    @test tapprox(test1, test3, atol=7e-13)
end
