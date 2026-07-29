# Combine or concatenate multiple formula objects

This small helper function allows to combine multiple
[`formula()`](https://rdrr.io/r/stats/formula.html) objects into one. In
the case of duplicate variable entries, only the unique ones are used.

## Usage

``` r
combine_formulas(..., combine = "both", env = parent.frame())
```

## Arguments

- ...:

  Any number [`formula`](https://rdrr.io/r/stats/formula.html) objects
  in "LHS ~ RHS" format, also supporting
  [`character`](https://rdrr.io/r/base/character.html) strings.

- combine:

  [`character`](https://rdrr.io/r/base/character.html) on whether LHS
  and RHS duplicates are to be removed. Can be set to either `"lhs"`,
  `"rhs"` or `"both"` (Default).

- env:

  A new environment of the formula `(def=parent.frame())`.

## Value

A formula as `cbind(lhs_1, lhs_2, ...) ~ rhs_1 + rhs_2 + ...` or
`lhs ~ rhs_1 + rhs_2` in case of identical LHS (see examples).

## Details

Use "y ~ 0" to specify a stand alone LHS.

## Note

This likely won't work for interaction terms (such as `*` or `:`).

## Examples

``` r
# Combine everything (default)
combine_formulas(observed ~ rainfall + temp, observed ~ rainfall + forest.cover)
#> observed ~ rainfall + temp
#> <environment: 0x559d0fd22e40>
# Combine only LHS
combine_formulas(observed ~ rainfall + temp, observed ~ rainfall + forest.cover, combine = "lhs")
#> observed ~ rainfall + temp + rainfall + forest.cover
#> <environment: 0x559d0fd22e40>
```
