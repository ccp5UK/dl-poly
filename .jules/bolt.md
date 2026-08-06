## 2024-05-24 - [Fortran Performance Optimization]
**Learning:** Float exponentiation (`x**2.0_wp`) is significantly slower than integer exponentiation (`x**2`) because float exponentiation invokes an expensive internal `exp(y * log(x))` library call. It also poses stability issues if negative bases are encountered.
**Action:** Consistently replace occurrences of floating point exponentiation like `**2.0_wp` or `**3.0_wp` with integer literals like `**2` or `**3` in computationally intensive sections.
