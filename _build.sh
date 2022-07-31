#!/bin/sh

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::bs4_book')"


git add --all
git commit -m "update book"
git push origin master
