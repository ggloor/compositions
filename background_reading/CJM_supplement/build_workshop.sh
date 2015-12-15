#!/bin/bash

set -e

FILE="workshop"
FIGS="figure"

rm -rf "${FIGS}"
Rscript -e "library(knitr); knit('${FILE}.Rnw')"
latexmk -xelatex -bibtex -pv "${FILE}"
latexmk -c "${FILE}"
