## 2024-05-18 - Avoid Floating Point Exponentiation in Fortran
**Learning:** Evaluating exponents using `**2.0_wp` invokes `exp(y * log(x))` from the math library. It is considerably slower than integer exponentiation (e.g., `**2`). Further, computing floating-point exponentiation for negative bases will cause NaN errors or mathematical exceptions. Thus, `x**2` or `x**3` is significantly faster and more robust than `x**2.0_wp` or `x**3.0_wp`.
**Action:** Consistently replace `**2.0_wp` with `**2` in performance-critical Fortran routines, including potential energy and integral evaluations.
