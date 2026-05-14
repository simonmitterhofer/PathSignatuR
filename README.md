# PathSignatuR

Truncated path signatures for multi-dimensional discrete paths. Built on
the exact piecewise-linear formulation (segment-wise tensor exponentials
and Chen products), with an Rcpp kernel and a pure-R reference for
verification.

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
s
#         1          2        1,1        1,2        2,1        2,2 
# 1.000  -0.83       1.42      0.34      -0.59      0.31       1.01

# Without the empty-word entry
signature(X, depth = 2, includeLevelZero = FALSE)

# Prepend a time channel (normalised to [0, 1])
signature(timeAugment(X), depth = 2)

# Inspect the word ordering
enumerateWords(dim = 2, depth = 2)
```

## Batch and path-valued variants

```r
# Many paths at once — precompute is shared across the batch
paths <- list(
  a = matrix(rnorm(40), 20, 2),
  b = matrix(rnorm(60), 30, 2),
  c = matrix(rnorm(50), 25, 2)
)
signatureBatch(paths, depth = 2)         # 3 x 7 matrix, row names a/b/c

# Running signature at every time step
sp <- signaturePath(X, depth = 2)        # T x 7 matrix
sp[1, ]                                  # identity (1, 0, ..., 0)
sp[nrow(sp), ]                           # same as signature(X, 2)
```

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

## Conventions

- **Path layout.** `X` is a `T x d` numeric matrix: rows are ordered
  observations, columns are channels. A numeric vector is treated as a
  1D path.
- **Word ordering.** Level-major, lexicographic within level (last
  letter varies fastest). The empty word sits at position 1.
- **Output.** `signature()` returns the terminal value of the truncated
  signature as a named numeric vector. `signaturePath()` returns the
  full running signature; `signatureBatch()` does many paths at once.
- **Time.** Not automatic. Use `timeAugment()` explicitly when you want
  a time channel.
- **Names.** Default `sep = ","` (`"1,1"`, `"1,2"`, ...). For
  `d <= 9` you may prefer `sep = ""` for compact `"11"`, `"12"`. For
  `d >= 10` keep the comma.

## Mathematical correctness

The Rcpp kernel uses segment-wise tensor exponentials and Chen products,
so closed forms are recovered to machine precision rather than
approximated:

```r
# Linear path: signature term at level k equals (b - a)^k / k!
x <- seq(0, 2, length.out = 100)
signature(x, depth = 4)
# 1.0000000 2.0000000 2.0000000 1.3333333 0.6666667
2^(0:4) / factorial(0:4)
# 1.0000000 2.0000000 2.0000000 1.3333333 0.6666667
```

## Roadmap

- v0.4: `tensorProduct()`, `chenProduct()`, `shuffleProduct()`,
  `tensorInverse()`, `visibilityTransform()`
- v0.5: `logSignature()`
- v1.0: stable API freeze

## License

MIT.
