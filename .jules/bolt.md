## 2024-05-24 - Integer vs Floating-Point Exponentiation in Fortran
**Learning:** Fortran evaluates floating-point exponents (e.g., `x**2.0_wp`) using expensive math library function calls like `exp(y * log(x))`. Furthermore, evaluating negative bases with floating-point exponents can cause NaN errors or exceptions in Fortran.
**Action:** Prefer integer exponentiation (e.g., `x**2` or `x**3`) over floating-point exponentiation for whole number powers to improve performance and code robustness.
