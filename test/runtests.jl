using SvdUpdate
using Test

#m=number of rows
#n=number of columns
#r=rank of matrix undergoing update
#p=rank of the update matrix

@testset "SvdUpdate.jl" begin
    #Write your tests here.
    #Test Case 1: Update to Square Rank Deficient Matrix
    m = 70; n = 70
    p = 10
    r = 5
    Y₁=rand(Complex{Float64},m,r)*rand(Complex{Float64},r,n)
    Y₂=rand(Complex{Float64},m,p)
    Y₃=rand(Complex{Float64},n,p)
    Ŷ=svd(Y₁)
    Ỹ=svdupdate(Ŷ,Y₂,Y₃)
    @test Ỹ.U*Diagonal(Ỹ.S)*Ỹ.V'≈Y₁+Y₂*Y₃'
    
    #Test Case 2: Update to Positive Definite Matrix
    m=70;n=70;r=70
    p=40
    Eᵩ=rand(m,n)
    E₁=Hermitian(Eᵩ)*Eᵩ
    E₂=rand(m,p)
    E₃=rand(n,p)
    Ê=svd(E₁)
    Ẽ=svdupdate(Ê,E₂,E₃)
    @test Ẽ.U*Diagonal(Ẽ.S)*Ẽ.V'≈E₁+E₂*E₃'
    
    #Test Case 3: Update to Full Rank Square Matrix
    m=80;n=80;r=80
    p=40
    Z₁=rand(Complex{Float64},m,n)
    Z₂=rand(Complex{Float64},m,p)
    Z₃=rand(Complex{Float64},n,p)
    Ẑ=svd(Z₁)
    Z̃=svdupdate(Ẑ,Z₂,Z₃)
    @test Z̃.U*Diagonal(Z̃.S)*Z̃.V'≈Z₁+Z₂*Z₃'
    
    #Test Case 4: Update to Hermitian Matrix
    m=65;n=65;r=40
    p=10
    T₁=Hermitian(rand(Complex{Float64},m,r)*rand(Complex{Float64},r,n))
    T₂=rand(m,p)
    T₃=rand(n,p)
    T̂=svd(T₁)
    T̃=svdupdate(T̂,T₂,T₃)
    @test T̃.U*Diagonal(T̃.S)*T̃.V'≈T₁+T₂*T₃'
    
    #Test Case 5: Update to Square Rank Deficient Matrix
    m=155;n=100;r=100
    p=1
    L₁=rand(m,r)*rand(r,n)
    L₂=rand(m,p)
    L₃=rand(n,p)
    L̂=svd(L₁)
    L̃=svdupdate(L̂,L₂,L₃)
    @test L̃.U*Diagonal(L̃.S)*L̃.V'≈L₁+L₂*L₃'
    
    #Test Case 6: Update to Tall Rank Deficient Matrix
    m=110;n=70;r=90
    p=40
    J₁=rand(Complex{Float64},m,r)*rand(Complex{Float64},r,n)
    J₂=rand(Complex{Float64},m,p)
    J₃=rand(Complex{Float64},n,p)
    Ĵ=svd(J₁)
    J̃=svdupdate(Ĵ,J₂,J₃)
    @test J̃.U*Diagonal(J̃.S)*J̃.V'≈J₁+J₂*J₃'    
    
    #Test Case 7: Update to Tall Full Rank Matrix
    m=170;n=100;r=100
    p=15
    W₁=rand(Complex{Float64},m,n)
    W₂=rand(Complex{Float64},m,p)
    W₃=rand(Complex{Float64},n,p)
    Ŵ=svd(W₁)
    W̃=svdupdate(Ŵ,W₂,W₃)
    @test W̃.U*Diagonal(W̃.S)*W̃.V'≈W₁+W₂*W₃'
    
    #Test Case 8: Update to Tall Full Rank Matrix with HIGH condition number
    m=50;n=45
    p=5
    G₁=rand(Complex{Float64},m,n)
    G₁[:,n]=G₁[:,1]+1e-10*ones(m)
    G₂=rand(Complex{Float64},m,p)
    G₃=rand(Complex{Float64},n,p)
    Ĝ=svd(G₁)
    G̃=svdupdate(Ĝ,G₂,G₃)
    G̃.U*Diagonal(G̃.S)*G̃.V'≈G₁+G₂*G₃'
    
    #Test Case 9: Update to Square matrix with HIGH condition number
    m=150;n=150;p=60
    H₁=rand(Complex{Float64},m,n)
    H₁[:,n]=H₁[:,1]+1e-10*ones(m)
    H₂=rand(Complex{Float64},m,p)
    H₃=rand(Complex{Float64},n,p)
    Ĥ=svd(H₁)
    H̃=svdupdate(Ĥ,H₂,H₃)
    H̃.U*Diagonal(H̃.S)*H̃.V'≈H₁+H₂*H₃'

end

#Creation of High Condition Number Matrices
C=rand(m,n)
C[:,n]=C[:,1]+1e-10*ones(m)
#this creates a matrix that is "Almost" rank deficient with a HUGE condition number
#make the final column "Almost" equal to the first column by redefining it to the first column plus a very tiny offset O(1e-10)










