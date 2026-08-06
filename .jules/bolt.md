## 2026-08-01 - [Integer Exponentiation Optimization]
**Learning:** In Fortran, using floating-point exponents (e.g., `x**2.0_wp` or `x**3.0_wp`) is a performance anti-pattern. It forces the compiler to use expensive math library calls like `exp(y * log(x))` instead of simple multiplications (`x * x`). This can also lead to NaN errors or exceptions with negative bases.
**Action:** Always prefer integer exponentiation (`x**2`, `x**3`) over floating point exponentiation when the exponent is a whole number.
