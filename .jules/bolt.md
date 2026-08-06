## 2024-06-17 - Floating Point Exponentiation in Fortran
**Learning:** In Fortran, floating-point exponentiation (e.g., `x ** 2.0_wp`) can be significantly slower than integer exponentiation (e.g., `x ** 2`) because the former often invokes the math library's `pow` or `exp(log())` functions, whereas the latter translates to direct multiplication.
**Action:** Always prefer integer exponents for powers like 2, 3, 4, etc., to avoid unnecessary overhead in hot paths like integrators and potential calculations.
