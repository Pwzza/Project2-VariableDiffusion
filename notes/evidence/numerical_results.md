# Generated numerical evidence

Run mode: full. See data/metadata.json.

## Smooth model: selected fits

| Method | Norm | Order | Fitted N | Max log deviation |
|---|---|---:|---|---:|
| fdm_arithmetic | l2 | 1.99989 | [32, 64, 128] | 2.94e-05 |
| fdm_arithmetic | nodal_linf | 1.99984 | [32, 64, 128] | 0.000135 |
| fem_fitted | l2 | 1.99898 | [32, 64, 128] | 0.000283 |
| fem_fitted | nodal_linf | 1.99939 | [32, 64, 128] | 7.23e-05 |
| fem_fitted | energy | 0.99956 | [32, 64, 128] | 0.000121 |

## Layered model: finest mesh

| Method | N | Nodal max error | Relative slice flux error | Boundary violation |
|---|---:|---:|---:|---:|
| fdm_arithmetic | 255 | 0.0023829 | 0.0031028 | 0 |
| fdm_resistance | 255 | 1.9995e-13 | 5.4746e-12 | 0 |
| fem_fitted | 255 | 5.4956e-13 | 1.1696e-11 | 0 |

## Acceptance gates

- smooth_l2_orders: PASS
- smooth_algebraic_separation: PASS
- quadrature_separation: PASS
- symmetric: PASS
- positive_eigenvalues: PASS
- layered_exact_nodes: PASS
- layered_exact_flux: PASS
- layered_bounds: PASS

Residual correction is a numerical algebraic-error estimate, not a certified bound.
Timings are single-run diagnostics; no controlled speed ranking is claimed.
