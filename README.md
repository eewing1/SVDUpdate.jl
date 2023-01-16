# SVDUpdate

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://eewing1.github.io/SVDUpdate.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://eewing1.github.io/SVDUpdate.jl/dev/)
[![Build Status](https://github.com/eewing1/SVDUpdate.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eewing1/SVDUpdate.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/eewing1/SVDUpdate.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/eewing1/SVDUpdate.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


This is a package based on the algorithm/identities created by Brand in the paper "Fast Low Rank Modifications of the Thin Singular Value Decomposition." The rough idea is as follows:

Given some matrix $X\in R^{mxn}$, computing the svd of the matrix is an expensive, an $O(n^3)$ operation for an n by n matrix. Suppose we already have the singular value decomposition $X=U\Sigma V^{H}$. Suppose we make some additive modification to the original matrix X, such as:
$$X+AB^{H}$$
Here, $A\in R^{mxq}$, anad $B\in{nxq}$
