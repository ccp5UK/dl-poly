#!/bin/bash
# make latexpdf generates a malformed dl_poly.ind file, resulting
# in a failed LaTeX build. Manually running pdflatex is a work around.
make latex
(cd build/latex; pdflatex dl_poly.tex && pdflatex dl_poly.tex)
