module SvdUpdate

using LinearAlgebra
# Write your package code here.
function svdupdate(F::SVD, A, B)
    U=F.U
    S=Diagonal(F.S)
    V=F.V

    σ₁=F.S[1]
    ϵₘ=eps()
    ê=σ₁*ϵₘ

    #remove all the left/right sing. vectors that have tiny sing. values
    #r=minimum(findall(svd(X).S.<10^-10
    r=findall(F.S.<ê)

    #test if matrix is rank deficient or full rank
    #r is an array of singular values that are TINY
    if isempty(r)
        r=size(F.U,2)
    else
        r=minimum(r)-1
    end

    Uᵣ=U[:,1:r]
    Vᵣ=V[:,1:r]
    Sᵣ=S[1:r,1:r]

    #compute the Q, R matrices
    Qₐ,Rₐ=qr(A-Uᵣ*Uᵣ'*A)
    Qᵦ,Rᵦ=qr(B-Vᵣ*Vᵣ'*B)

    #convert "full" Q matrices into "thin" matrices that match dimensions of A,B
    #the Q that is returned by the qr function is not what we're looking for
    Qₐ=Matrix(Qₐ)
    Qᵦ=Matrix(Qᵦ)

    #build the K matrix
    K₁=Sᵣ+Uᵣ'*A*B'*Vᵣ
    K₂=Uᵣ'*A*Rᵦ'
    K₃=Rₐ*B'*Vᵣ
    K₄=Rₐ*Rᵦ'
    K=[K₁ K₂;K₃ K₄]

    #build arrays whose entries are diagonals of Ra,Rb
    Rₐ₁=abs.(diag(Rₐ))
    Rₐ₂=abs.(diag(Rᵦ))

    #size of the array Ra, Rb
    rₐ₁=size((findall(Rₐ₁.<ê)),1)
    rₐ₂=size((findall(Rₐ₂.<ê)),1)

    #find dimensions of K, truncate K
    mₖ,nₖ=size(K)
    K=K[1:(mₖ-rₐ₁),1:(nₖ-rₐ₂)]

    #number of columns of Qa, Qb
    qₐ=size(Matrix(Qₐ),2)
    qᵦ=size(Matrix(Qᵦ),2)


    if qₐ==rₐ₁
        Qₐ=Matrix{Float64}(undef,size(U,1),0)
    elseif size(Qₐ,2)>rₐ₁
        Qₐ=Qₐ[:,1:(qₐ-rₐ₁)]
    end

    if qᵦ==rₐ₂
        Qᵦ=Matrix{Float64}(undef,size(V,1),0)
    elseif size(Qᵦ,2)>rₐ₂
        Qᵦ=Qᵦ[:,1:(qᵦ-rₐ₂)]
    end



    #compute the SVD of K
    Uₖ,Sₖ,Vₖ=svd(K)

    SVD(([Uᵣ Qₐ]*Uₖ), Sₖ, ([Vᵣ Qᵦ]*Vₖ)');
end
end # module
