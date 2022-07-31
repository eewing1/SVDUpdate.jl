using SvdUpdate
using Test

@testset "SvdUpdate.jl" begin
    #Write your tests here.
    #Test Case 1: Update to Square Rank Deficient Matrix
    m = 70; n = 70
    p = 10
    r = 5
    Y₁=rand(Complex{Float64},150,120)*rand(Complex{Float64},120,150)
    Y₂=rand(Complex{Float64},150,30)
    Y₃=rand(Complex{Float64},150,30)
    Ŷ=svd(Y₁)
    Ỹ=totalbrand6(Ŷ,Y₂,Y₃)
    @test Ỹ.U*Diagonal(Ỹ.S)*Ỹ.V'≈Y₁+Y₂*Y₃'
    
    #Test Case 2: Update to Positive Definite Matrix
    m=70;n=70;r=70
    p=40
    Eᵩ=rand(m,n)
    E₁=Hermitian(Eᵩ)*Eᵩ
    E₂=rand(m,p)
    E₃=rand(n,p)
    Ê=svd(E₁)
    Ẽ=totalbrand6(Ê,E₂,E₃)
    @test Ẽ.U*Diagonal(Ẽ.S)*Ẽ.V'≈E₁+E₂*E₃'
    
    #Test Case 3: Update to Full Rank Square Matrix
    m=80;n=80;r=80
    p=40
    Z₁=rand(Complex{Float64},m,n)
    Z₂=rand(Complex{Float64},m,p)
    Z₃=rand(Complex{Float64},n,p)
    Ẑ=svd(Z₁)
    Z̃=totalbrand6(Ẑ,Z₂,Z₃)
    @test Z̃.U*Diagonal(Z̃.S)*Z̃.V'≈Z₁+Z₂*Z₃'
    
    #Test Case 4: Update to Hermitian Matrix
    m=65;n=65;r=40
    p=10
    T₁=Hermitian(rand(Complex{Float64},m,r)*rand(Complex{Float64},r,n))
    T₂=rand(m,p)
    T₃=rand(n,p)
    T̂=svd(T₁)
    T̃=totalbrand6(T̂,T₂,T₃)
    @test T̃.U*Diagonal(T̃.S)*T̃.V'≈T₁+T₂*T₃'
    
    #Test Case 5: Update to Square Rank Deficient Matrix
    m=155;n=100;r=100
    p=1
    L₁=rand(m,r)*rand(r,n)
    L₂=rand(m,p)
    L₃=rand(n,p)
    L̂=svd(L₁)
    L̃=totalbrand6(L̂,L₂,L₃)
    @test L̃.U*Diagonal(L̃.S)*L̃.V'≈L₁+L₂*L₃'
    
    #Test Case 6: Update to Tall Rank Deficient Matrix
    m=110;n=70;r=90
    p=40
    J₁=rand(Complex{Float64},m,r)*rand(Complex{Float64},r,n)
    J₂=rand(Complex{Float64},m,p)
    J₃=rand(Complex{Float64},n,p)
    Ĵ=svd(J₁)
    J̃=totalbrand6(Ĵ,J₂,J₃)
    @test J̃.U*Diagonal(J̃.S)*J̃.V'≈J₁+J₂*J₃'    

end
