## YYYY-MM-DD - [Title]

## 2026-06-18 - Prefer integer exponentiation over floating-point
**Learning:** In Fortran, evaluating mathematical expressions with small integer powers (like x**2 or x**3) using floating point representations of the exponent (like `x**2.0_wp`) is a common anti-pattern in this codebase. Floating-point exponentiation is evaluated using expensive math library calls like `exp(y * log(x))` which slows down computationally intensive inner loops. Moreover, evaluating negative bases with floating-point exponents can cause NaN errors or exceptions in Fortran.
**Action:** Always prefer integer exponentiation (`x**2`) instead of floating-point exponents (`x**2.0`) for small integer powers to ensure they are evaluated efficiently via simple multiplications.
