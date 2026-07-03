## 2024-05-14 - Integer Exponentiation over Floating Point
**Learning:** In Fortran, integer exponentiation (e.g. `x**2`) is evaluated differently from floating-point exponentiation (e.g. `x**2.0_wp`). Floating point powers involve calling the math library (`exp(y * log(x))`), which can be computationally expensive and unstable for negative bases.
**Action:** When a power is known to be a small integer, always write it as an integer literal rather than a real variable to avoid performance penalty and ensure robust mathematical behavior.
