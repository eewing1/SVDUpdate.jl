# SvdUpdate

[![Build Status](https://github.com/eewing1/SvdUpdate.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/eewing1/SvdUpdate.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Here we seek to build a Julia package that performs low rank updates to a matrix, given its singular value decomposition. Our inspiration is given by an algorithm in "Fast Methods of for QR Decomposition" by Matthew Brand, linked below.

https://www.sciencedirect.com/science/article/pii/S0024379505003812


## Additive Modications to the SVD

let $X\in R^{mxn}$, rank r, thin SVD $X=U\Sigma{V^H}$ <br>
let $A\in R^{mxc}$, $B\in R^{nxc}$ <br>
What is the SVD of below sum, which is a rank "c" update to matrix X? <br>
Can we express the SVD of this sum as updates to $U,S,V$?

$$
X+A{B^H}=
\begin{pmatrix} U&A\end{pmatrix}
\begin{pmatrix} \Sigma & 0\\ 0&I\end{pmatrix}
\begin{pmatrix} V&B\end{pmatrix}^H
$$

###  (i) General rank case
#### a) Matrices pertaining to output space
- $A-U{U^H}{A}$ is a matrix that gives us the part of $C(A)$ that is orthogonal to $C(U)$
        - what information does A have that is not already contained in $C(U)$

- Perform a <font color=red>**QR Decomposition**</font> on this matrix $A-U{U^H}A$ 
    - Let $Q_a$ be an <font color=blue>orthogonal basis</font> of the column space of $A-U{U^H}{A}$
    - let $R_a$ be an upper triangular matrix obtained via this QR decomposition

The relationship between these matrices $U$, $A$, $Q_a$, $R_a$, etc., is summarized below:
$
\begin{pmatrix}U&A\end{pmatrix}
=
\begin{pmatrix}U&Q_a\end{pmatrix}
\begin{pmatrix}I&{U^H}{A}\\0& R_a\end{pmatrix}
$

#### b) Matrices pertaining to input space
- $B-V{V^H}B$ is a matrix that gives the part of $C(B)$ that is orthogonal to $C(V)$
    - Perform a <font color=red>**QR Decomposition**</font> on this matrix $B-V{V^H}B$ 

The relationship between these matrices $B$, $V$, $Q_b$, $R_b$ is summarized as follows:

$
\begin{pmatrix}V&B\end{pmatrix}
=
\begin{pmatrix}V&Q_b\end{pmatrix}
\begin{pmatrix}I&{V^H}{B}\\0& R_b\end{pmatrix}
$

Remember: our earlier equation was:
$$
X+{A}{B^H}=
\begin{pmatrix} U&A\end{pmatrix}
\begin{pmatrix} \Sigma & 0\\0&I\end{pmatrix}
\begin{pmatrix} V&B\end{pmatrix}^H
$$
