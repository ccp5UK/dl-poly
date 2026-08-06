## 2024-06-20 - Fortran Exponentiation Optimization
**Learning:** In Fortran codebases like DL_POLY_4, using floating-point exponentiation (e.g., `x**2.0_wp`) forces the compiler to use expensive math functions like `exp(y * log(x))`. Integer exponentiation (`x**2`) compiles to simple multiplication, which is significantly faster without losing precision.
**Action:** When reviewing Fortran calculations inside hot loops (like integration or potentials), actively hunt for and replace `**N.0_wp` with `**N`.
