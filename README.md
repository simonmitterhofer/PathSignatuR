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

## Conventions

- **Path layout.** `X` is a `T x d` numeric matrix: rows are ordered
  observations, columns are channels. A numeric vector is treated as a
  1D path. v0.1 is single-path only.
- **Word ordering.** Level-major, lexicographic within level (last
  letter varies fastest). The empty word sits at position 1.
- **Output.** `signature()` returns the **terminal** value of the
  truncated signature as a named numeric vector. The full path-valued
  signature is computed internally; an exported `signaturePath()` is
  on the v0.2 roadmap.
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

- v0.2: `signatureBatch()`, `signaturePath()`
- v0.3: `leadLag()`, `basepointAugment()`, `visibilityTransform()`,
  `pathInterpolate()`
- v0.4: `tensorProduct()`, `chenProduct()`, `shuffleProduct()`
- v0.5: `logSignature()`

## License

MIT.
