#!/bin/sh

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"


git add --all
git commit -m "update book"
git push origin master
