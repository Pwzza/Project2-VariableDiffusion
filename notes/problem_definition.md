# Problem definition: Direction 3

Status: proposed conventions and derived predictions, written before project solves.

## Part 1: continuous problem and purpose

Model steady conduction through a heterogeneous material. Let temperature be `u`, conductivity be `a`, and heat flux be `q = -a grad(u)`. Balance gives `div(q) = f`, hence

\[
-\nabla\cdot(a\nabla u)=f\quad\text{in }\Omega=(0,1)^2,
\qquad u=g\quad\text{on }\partial\Omega,
\qquad 0<a_{\min}\le a\le a_{\max}<\infty.
\]

The conductivity is scalar and isotropic. We study a steady, linear, prescribed-temperature problem. No transient, nonlinear, anisotropic, contact-resistance, or experimental-material claims are implied.

At a perfect interface without an interface source, `[u]=0` and `[q dot n]=0` using the **same oriented normal** on both sides. In outward normals of the two material subdomains the two fluxes instead sum to zero. A jump in the temperature gradient is physically expected.

**Method choice: conservative Cartesian FDM.** The square and vertical interface make shared face fluxes explicit, allow transparent boundary elimination, and isolate arithmetic versus resistance averaging. Use sparse direct `spsolve` for the baseline and residual-stopped Gauss-Seidel for the selected depth comparison. Add CG/AMG only after that comparison is complete.

Two required benchmarks:

1. Smooth: `a=1+x`, `u*=sin(pi*x)sin(pi*y)`, `g=0`, forcing derived in the experiment plan.
2. Layered: `a=1` left of `x=1/2`, `a=10` right, `f=0`, `g=u*` on **all four sides**, with the continuous piecewise-linear exact profile from the pack. Never replace the top/bottom data by zero or by insulation in the core test.

## Part 2: discrete operator and conventions

`N` always means **number of intervals per axis**, `h=1/N`, and `n=N-1` interior nodes per axis. Node `(i,j)` has `(x,y)=(i*h,j*h)`, `i,j=1,...,N-1`. Interior vector index is

\[
k=(i-1)n+(j-1).
\]

Use `meshgrid(..., indexing='ij')`, C-order flatten/reshape; axis 0 is x. The course setup checker uses `n` for interior nodes, so its `n=15` corresponds to project `N=16`.

Store each horizontal/vertical edge conductivity once, then use it in both adjoining rows:

\[
(A_hU)_{ij}=h^{-2}\{a_E(U_{ij}-U_{i+1,j})+
a_W(U_{ij}-U_{i-1,j})+a_N(U_{ij}-U_{i,j+1})+
a_S(U_{ij}-U_{i,j-1})\}.
\]

Diagonal is the sum of all four positive face values divided by `h^2`. An interior neighbor contributes `-a_face/h^2`; a Dirichlet neighbor contributes `+a_face*g/h^2` to `b`, while its diagonal contribution remains. Unknowns are interior nodes only. Assemble directly into COO and convert once to CSR. Use `A @ U`, never an ambiguous `A * U`.

For smooth data use the arithmetic mean of endpoint conductivities. Across an interface at known distances `dL+dR=h`, use

\[
a_{\mathrm{face}}=\frac{h}{d_L/a_L+d_R/a_R}.
\]

At a bisected edge this is `2*aL*aR/(aL+aR)`. Arithmetic averaging at that jump is the intentional adverse comparison only. For core odd-interval meshes the interface never coincides with a node, so no arbitrary value of `a(1/2)` is needed.

For a zero-boundary vector v extended by zero to boundary nodes,

\[
h^2v^TA_hv=\sum_{\{p,q\}\in\mathcal E}a_{pq}(v_p-v_q)^2>0\quad(v\ne0),
\]

where every undirected interior-interior or interior-boundary edge is counted once. This proves SPD: positivity plus connectivity to the Dirichlet boundary removes the constant nullspace. Shared weights prove symmetry. The nonpositive off-diagonals and irreducible weak diagonal dominance, strict next to the boundary, give an M-matrix. Thus with `f=0` the solution lies between boundary extrema. Positive contrast preserves these properties.

This is a proof about the operator. Numerical asymmetry and smallest eigenvalues corroborate it; random Rayleigh quotients alone do not prove it. Replacing boundary rows by identity while retaining their columns is excluded because it destroys symmetry.

## Hand-check before implementation

For `N=3`, constant `a=1`, four interior unknowns in declared order:

\[
A=9\begin{pmatrix}4&-1&-1&0\\-1&4&0&-1\\-1&0&4&-1\\0&-1&-1&4\end{pmatrix}.
\]

For row `(i,j)=(1,1)` under nonzero boundary values,
`b[0]=f(h,h)+9*g(0,h)+9*g(h,0)`. For variable `a`, replace each `9` by that edge's `a_face/h^2` and the diagonal by the sum of its four faces. Check both x-neighbor stride `n` and y-neighbor stride `1` without allowing a wrap between rows.

For the bisected interface with `aL=1,aR=10`, the crossing face is `20/11`; the intentional arithmetic variant is `11/2`. Both rows adjacent to this face must use the identical stored value.
