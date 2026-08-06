## 2024-05-15 - Fortran Floating Point Exponentiation Penalty
**Learning:** In Fortran, raising a number to a real/floating-point power (e.g. `x**2.0_wp`) is significantly slower than integer exponentiation (`x**2`), as the compiler often uses `exp(y * log(x))` behind the scenes. This causes a massive performance hit when repeated in inner loops (like force or energy calculations).
**Action:** Always scan for and replace instances like `**2.0` or `**2.0_wp` with `**2` in computationally heavy `.F90` routines.
