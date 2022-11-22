"""
    svdupdate(F, A, B)

Update the SVD.
"""
function svdupdate(F::SVD, A, B)
    U = F.U
    S = Diagonal(F.S)
    V = F.V

    σ₁ = F.S[1]
    ϵₘ = eps()
    ê = σ₁ * ϵₘ

    #remove all the left/right sing. vectors that have tiny sing. values
    #count function counts number of elements of F.S that satisfy condition <ehat
    r = count(<(ê), F.S)

    #test if matrix is rank deficient or full rank
    #r is an array telling us how many rank deficient singular values
    #how many singular vectors should be chopped off U, V

    uᵣ = size(U, 2) - r
    vᵣ = size(V, 2) - r
    sᵣ = size(S, 2) - r
    #Vᵣ=V[:,1:size(V,2)-r]
    # Sᵣ=S[1:s,1:s]
    Uᵣ = @views U[:, 1:uᵣ]
    Vᵣ = @views V[:, 1:vᵣ]
    Sᵣ = @views S[1:sᵣ, 1:sᵣ]

    #compute the Q, R matrices
    Qₐ, Rₐ = qr(A - (Uᵣ) * (Uᵣ)' * A)
    Qᵦ, Rᵦ = qr(B - (Vᵣ) * (Vᵣ)' * B)

    #convert "full" Q matrices into "thin" matrices that match dimensions of A,B
    #the Q that is returned by the qr function is not what we're looking for
    qₐ = size(Rₐ, 2)
    qᵦ = size(Rᵦ, 2)


    #build the K matrix
    K₁ = (Sᵣ) + (Uᵣ)' * A * B' * (Vᵣ)
    K₂ = (Uᵣ)' * A * Rᵦ'
    K₃ = Rₐ * B' * (Vᵣ)
    K₄ = Rₐ * Rᵦ'
    K = [K₁ K₂; K₃ K₄]

    #build arrays whose entries are diagonals of Ra,Rb
    rₐ₁ = count(<(ê), abs(Rₐ[i, i]) for i in minimum(axes(Rₐ)))
    rₐ₂ = count(<(ê), abs(Rᵦ[i, i]) for i in minimum(axes(Rᵦ)))



    #find dimensions of K, truncate K
    mₖ, nₖ = size(K)

    kᵣ = (mₖ - rₐ₁)
    kₙ = (nₖ - rₐ₂)
    Kᵣ = (@views K[1:kᵣ, 1:kₙ])

    #number of columns of Qa, Q

    #compute the SVD of K
    Uₖ, Sₖ, Vₖ = svd(Kᵣ)
    SVD([Uᵣ (@views Qₐ[:, 1:qₐ-rₐ₁])] * Uₖ, Sₖ, ([Vᵣ (@views Qᵦ[:, 1:qᵦ-rₐ₂])] * Vₖ)')
end
