module SvdUpdate

function svdupdate(F::SVD, A, B)
    U=F.U
    S=Diagonal(F.S)
    V=F.V

    σ₁=F.S[1]
    ϵₘ=eps()
    ê=σ₁*ϵₘ

    #remove all the left/right sing. vectors that have tiny sing. values
    #count function counts number of elements of F.S that satisfy condition <ehat
    r=count(<(ê), F.S)

    #test if matrix is rank deficient or full rank
    #r is an array telling us how many rank deficient singular values
        #how many singular vectors should be chopped off U, V
    
    uᵣ= size(U,2)-r
    vᵣ=size(V,2)-r
    sᵣ=size(S,2)-r

    #compute the Q, R matrices
    Qₐ,Rₐ=qr(A-(@views U[:, 1:uᵣ])*(@views U[:, 1:uᵣ])'*A)
    Qᵦ,Rᵦ=qr(B-(@views V[:,1:vᵣ])*(@views V[:,1:vᵣ])'*B)

    #convert "full" Q matrices into "thin" matrices that match dimensions of A,B
    #the Q that is returned by the qr function is not what we're looking for
    Qₐ=Matrix(Qₐ)
    Qᵦ=Matrix(Qᵦ)
    

    #build the K matrix
    K₁=(@views S[1:sᵣ,1:sᵣ])+(@views U[:,1:uᵣ])'*A*B'*(@views V[:,1:vᵣ])
    K₂=(@views U[:,1:uᵣ])'*A*Rᵦ'
    K₃=Rₐ*B'*(@views V[:,1:vᵣ])
    K₄=Rₐ*Rᵦ'
    K=[K₁ K₂;K₃ K₄]

    #build arrays whose entries are diagonals of Ra,Rb
    rₐ₁=count(<(ê), abs(Rₐ[i,i]) for i in minimum(axes(Rₐ)))
    rₐ₂=count(<(ê), abs(Rᵦ[i,i]) for i in minimum(axes(Rᵦ)))

    

    #find dimensions of K, truncate K
    mₖ,nₖ=size(K)
    
    kᵣ=(mₖ-rₐ₁)
    kₙ=(nₖ-rₐ₂)
    @views K[1:kᵣ,1:kₙ]

    #number of columns of Qa, Qb
    qₐ=size(Matrix(Qₐ),2)
    qᵦ=size(Matrix(Qᵦ),2)


    if qₐ==rₐ₁
        Qₐ=Matrix{Float64}(undef,size(U,1),0)
    elseif qₐ>rₐ₁
        Qₐ=Qₐ[:,1:(qₐ-rₐ₁)]
    end

    if qᵦ==rₐ₂
        Qᵦ=Matrix{Float64}(undef,size(V,1),0)
    elseif qᵦ>rₐ₂
        Qᵦ=Qᵦ[:,1:(qᵦ-rₐ₂)]
    end



    #compute the SVD of K
    Uₖ,Sₖ,Vₖ=svd(@views K[1:kᵣ,1:kₙ])

   SVD([(@views U[:, 1:uᵣ]) Qₐ]*Uₖ, Sₖ, ([(@views V[:, 1:vᵣ]) Qᵦ]*Vₖ)')
    
end
end # module
