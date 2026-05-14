# PathSignatuR

Truncated path signatures for multi-dimensional discrete paths. Exact
piecewise-linear formulation (segment-wise tensor exponentials and Chen
products) with an Rcpp kernel, cached `(d, depth)` workspaces, and a
pure-R reference for verification.

## Install

```r
# from a local clone
devtools::install("path/to/PathSignatuR")

# or, once on GitHub
# remotes::install_github("yourname/PathSignatuR")
```

Requires a working C++ toolchain (Rtools on Windows, Xcode CLT on macOS).

## Quickstart

```r
library(PathSignatuR)

# 2D random walk
X <- matrix(cumsum(rnorm(200)), 100, 2)

# Truncated signature up to level 2 (named numeric vector)
s <- signature(X, depth = 2)

# Without the empty-word entry
signature(X, depth = 2, includeLevelZero = FALSE)

# Prepend a time channel (normalised to [0, 1])
signature(timeAugment(X), depth = 2)

# Inspect the word ordering
enumerateWords(dim = 2, depth = 2)
```

## Signatures and log-signatures

Both come in terminal, path-valued, and batch variants. All share the
cached `(d, depth)` workspace, so repeated calls at the same truncation
amortize the precompute.

```r
# Terminal signature / log-signature
signature(X, depth = 3)
logSignature(X, depth = 3)

# Running variants: T x p matrix, row t = signature of X[1:t, ]
signaturePath(X, depth = 3)
logSignaturePath(X, depth = 3)

# Batch over many paths
paths <- list(
  a = matrix(rnorm(40), 20, 2),
  b = matrix(rnorm(60), 30, 2),
  c = matrix(rnorm(50), 25, 2)
)
signatureBatch(paths, depth = 2)         # nPaths x p, rows named a/b/c
```

The kernel is C++ throughout. At `(d, N, T) = (3, 4, 100)`:
`signature` ≈ 2 ms, `signaturePath` ≈ 2 ms, `logSignature` ≈ 2 ms,
`logSignaturePath` ≈ 4 ms. `logSignaturePath` outpaces the naive
R loop over `logSignature` by ~130x.

## Lyndon-basis log-signatures

The log-signature lives in the free Lie algebra of rank `d` -- a strict
subspace of the truncated tensor algebra. Carrying all `sum(d^(0:N))`
tensor coordinates is redundant: Jacobi identity and bracket
antisymmetry impose linear relations. Project onto the canonical Lyndon
basis to get independent coordinates of Witt-formula dimension:
`(d, N) = (3, 4)` is **32** vs 121.

```r
# Independent Lyndon coordinates
beta <- logSignatureLyndon(X, depth = 4)
length(beta)                     # = length(lyndonWords(d = 2, 4))
names(beta)                      # "1" "2" "1,2" "1,1,2" "1,2,2" ...

# Running variant: T x q matrix
logSignaturePathLyndon(X, depth = 4)

# Inspect the Lyndon-word ordering
lyndonWords(dim = 2, depth = 4)

# Inverse: Lyndon coords back to tensor basis
expandLyndon(beta, dim = 2, depth = 4)
```

This is the *right* dimension reduction for log-signatures: algebraic,
not statistical. Apply on top before PCA, JL, or any other downstream
reduction.

## Path transforms

All transforms are path-in / path-out and compose cleanly with
`signature()` and friends.

```r
# Time augmentation: prepend a time channel (default unit-scaled to [0, 1])
timeAugment(X)

# Basepoint augmentation: prepend a reference point (default origin) to
# break translation invariance. Usually unnecessary for log-return work.
basepointAugment(X)

# Lead-lag (FHL): doubles channels into lead/lag copies. Level-2 cross-
# term S^{lead_i, lag_i} - S^{lag_i, lead_i} recovers the quadratic
# variation of channel i.
leadLag(X)

# Path interpolation: resample to a new grid (linear by default; "step"
# for left-continuous piecewise-constant). No splines by design --
# would break piecewise-linearity.
pathInterpolate(X, n = 50)
pathInterpolate(X, at = c(0, 0.25, 0.5, 0.75, 1))
```

## Tensor algebra primitives

Operations on truncated tensor algebra elements (vectors of length
`sum(d^(0:depth))` in `enumerateWords` order). All inputs and outputs
are named numeric vectors. The C++ kernel caches `(d, depth)`-dependent
workspaces.

```r
# Signatures concatenate via the tensor / Chen product
s1 <- signature(X[1:15, ], depth = 3)
s2 <- signature(X[15:30, ], depth = 3)
tensorProduct(s1, s2, depth = 3)         # == signature(X, depth = 3)
chenProduct(s1, s2, depth = 3)           # named alias for clarity

# Reverse path: tensor inverse
sRev <- tensorInverse(s1, depth = 3)     # == signature(X1 reversed)

# Shuffle product encodes the algebraic relations between signature
# coordinates: <S, u> * <S, v> = <S, u ⧢ v>
shuffleProduct(s1, s2, depth = 3)
```

## Conventions

- **Path layout.** `X` is a `T x d` numeric matrix: rows are ordered
  observations, columns are channels. A numeric vector is treated as a
  1D path.
- **Word ordering.** Level-major, lexicographic within level (last
  letter varies fastest). The empty word sits at position 1.
- **Output.** Terminal functions return named numeric vectors. Path
  variants return `T x p` numeric matrices. Batch variants return
  `nPaths x p`.
- **Time.** Not automatic. Use `timeAugment()` explicitly.
- **Names.** Default `sep = ","` (`"1,1"`, `"1,2"`, ...). For
  `d <= 9` you may prefer `sep = ""` for compact `"11"`, `"12"`. For
  `d >= 10` keep the comma.

## Mathematical correctness

The C++ kernel uses exact segment-wise tensor exponentials, so closed
forms are recovered to machine precision:

```r
# Linear path: signature term at level k equals (b - a)^k / k!
x <- seq(0, 2, length.out = 100)
signature(x, depth = 4)
# 1.0000000 2.0000000 2.0000000 1.3333333 0.6666667
2^(0:4) / factorial(0:4)
# 1.0000000 2.0000000 2.0000000 1.3333333 0.6666667
```

Verified against a pure-R reference (`.signatureRef`) in the test
suite (1300+ tests).

## Roadmap

- v0.7+: signature kernels (Salvi et al.), invertSignature, vignettes
- v1.0: stable API freeze

## License

MIT.
