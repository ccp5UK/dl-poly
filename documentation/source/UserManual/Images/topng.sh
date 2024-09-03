for i in *.svg; do inkscape $i --export-width=512 -o "${i%.*}.png"; done
