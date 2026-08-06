## 2024-06-25 - Avoid Floating Point Exponentiation
**Learning:** In this Fortran codebase, using floating-point exponents like `**2.0_wp` instead of integer exponents like `**2` forces the compiler to use expensive math library calls like `exp(y * log(x))` and makes it susceptible to `NaN` errors on negative bases.
**Action:** Always replace floating-point exponentiation with integer exponentiation for integers constants to get better performance.
