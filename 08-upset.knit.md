# UpSet plot {#upset-plot}



[UpSet plot](https://caleydo.org/tools/upset/) provides an efficient way to
visualize intersections of multiple sets compared to the traditional
approaches, i.e. the Venn Diagram. It is implemented in the [UpSetR
package](https://github.com/hms-dbmi/UpSetR) in R. Here we re-implemented
UpSet plots with the **ComplexHeatmap** package with some improvements.

## Input data {#input-data}

To represent multiple sets, the variable can be represented as:

1. A list of sets where each set is a vector, e.g.:


```r
list(set1 = c("a", "b", "c"),
     set2 = c("b", "c", "d", "e"),
     ...)
```

2. A binary matrix/data frame where rows are elements and columns are sets,
   e.g.:

```
  set1 set2 set3
h    1    1    1
t    1    0    1
j    1    0    0
u    1    0    1
w    1    0    0
...
```

E.g., for row `t`, it means, `t` is in set **set1**, not in set **set2**, and
in set **set3**. Note the matrix is also valid if it is a logical matrix.

If the variable is a data frame, the binary columns (only contain 0 and 1) and
the logical columns are only used.

Both formats can be used for making UpSet plots, users can still use
`list_to_matrix()` to convert from list to the binary matrix.


```r
lt = list(set1 = c("a", "b", "c"),
          set2 = c("b", "c", "d", "e"))
list_to_matrix(lt)
```

```
##   set1 set2
## a    1    0
## b    1    1
## c    1    1
## d    0    1
## e    0    1
```

You can also set the universal set in `list_to_matrix()`:


```r
list_to_matrix(lt, universal = letters[1:10])
```

```
##   set1 set2
## a    1    0
## b    1    1
## c    1    1
## d    0    1
## e    0    1
## f    0    0
## g    0    0
## h    0    0
## i    0    0
## j    0    0
```

3. The set can be genomic intervals, then it can only be represented as a list
   of `GRanges`/`IRanges` objects.


```r
list(set1 = GRanges(...),
	 set2 = GRanges(...),
	 ...)
```

## Mode {#upset-mode}

E.g. for three sets (**A**, **B**, **C**), all combinations of selecting
elements in the set or not in the set are as following:

```
A B C
1 1 1
1 1 0
1 0 1
0 1 1
1 0 0
0 1 0
0 0 1
```

A value of 1 means to select that set and 0 means not to select that set.
E.g., "1 1 0" means to select set A, B while not set C. Note there is no "0 0
0", because the background set is not of interest here. In following part of
this section, we refer **A**, **B** and **C** as **sets** and each combination
as **combination set**. The whole binary matrix is called **combination
matrix**.

The UpSet plot visualizes the **size** of each combination set. With the
binary code of each combination set, next we need to define how to calculate
the size of that combination set. There are three modes:

1. `distinct` mode: 1 means in that set and 0 means not in that set, then `1 1
   0` means a set of elements both in set **A** and **B**, while not in **C**
   (`setdiff(intersect(A, B), C)`). Under this mode, the seven combination
   sets are the seven partitions in the Venn diagram and they are mutually
   exclusive.

2. `intersect` mode: 1 means in that set and 0 is not taken into account,
   then, `1 1 0` means a set of elements in set **A** and **B**, and they can
   also in **C** or not in **C** (`intersect(A, B)`). Under this mode, the
   seven combination sets can overlap.

3. `union mode`: 1 means in that set and 0 is not taken into account. When
   there are multiple 1, the relationship is __OR__. Then, `1 1 0` means a set
   of elements in set **A** or **B**, and they can also in **C** or not in
   **C** (`union(A, B)`). Under this mode, the seven combination sets can
   overlap.

The three modes are illustrated in following figure:

<img src="08-upset_files/figure-html/unnamed-chunk-6-1.png" width="460.8" style="display: block; margin: auto;" />

## Make the combination matrix {#make-the-combination-matrix}

The function `make_comb_mat()` generates the combination matrix as well as
calculates the size of the sets and the combination sets. The input can be one
single variable or name-value pairs:


```r
set.seed(123)
lt = list(a = sample(letters, 5),
	      b = sample(letters, 10),
	      c = sample(letters, 15))
m1 = make_comb_mat(lt)
m1
```

```
## A combination matrix with 3 sets and 7 combinations.
##   ranges of combination set size: c(1, 8).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
##   x      100    1
##     x    010    3
##       x  001    8
## 
## Sets are:
##   set size
##     a    5
##     b   10
##     c   15
```

```r
m2 = make_comb_mat(a = lt$a, b = lt$b, c = lt$c)
m3 = make_comb_mat(list_to_matrix(lt))
```

`m1`, `m2` and `m3` are identical.

The mode is controlled by the `mode` argument:


```r
m1 = make_comb_mat(lt) # the default mode is `distinct`
m2 = make_comb_mat(lt, mode = "intersect")
m3 = make_comb_mat(lt, mode = "union")
```

The UpSet plots under different modes will be demonstrated in later sections.

When there are too many sets, the sets can be pre-filtered by the set sizes.
The `min_set_size` and `top_n_sets` are for this purpose. `min_set_size`
controls the minimal size for the sets and `top_n_sets` controls the number of
top sets with largest sizes.


```r
m1 = make_comb_mat(lt, min_set_size = 4)
m2 = make_comb_mat(lt, top_n_sets = 2)
```

The subsetting of the sets affects the calculation of the sizes of the
combination sets, that is why it needs to be controlled at the combination
matrix generation step. The subsetting of combination sets can be directly
performed by subsetting the matrix:


```r
m = make_comb_mat(lt)
m[1:4]
```

```
## A combination matrix with 3 sets and 4 combinations.
##   ranges of combination set size: c(1, 4).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
## 
## Sets are:
##   set size
##     a    5
##     b   10
##     c   15
```

`make_comb_mat()` also allows to specify the universal set so that the complement
set which contains elements not belonging to any set is also considered.


```r
m = make_comb_mat(lt, universal_set = letters)
m
```

```
## A combination matrix with 3 sets and 8 combinations.
##   ranges of combination set size: c(1, 8).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
##   x      100    1
##     x    010    3
##       x  001    8
##          000    6
## 
## Sets are:
##          set size
##            a    5
##            b   10
##            c   15
##   complement    6
```

The universal set can be smaller than the union of all sets, then for each
set, only the intersection to universal set is considered.


```r
m = make_comb_mat(lt, universal_set = letters[1:10])
m
```

```
## A combination matrix with 3 sets and 5 combinations.
##   ranges of combination set size: c(1, 3).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x    110    1
##   x   x  101    1
##     x x  011    2
##       x  001    3
##          000    3
## 
## Sets are:
##          set size
##            a    2
##            b    3
##            c    6
##   complement    3
```

If you already know the size of the complement size, you can directly set 
`complement_size` argument.


```r
m = make_comb_mat(lt, complement_size = 5)
m
```

```
## A combination matrix with 3 sets and 8 combinations.
##   ranges of combination set size: c(1, 8).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
##   x      100    1
##     x    010    3
##       x  001    8
##          000    5
## 
## Sets are:
##          set size
##            a    5
##            b   10
##            c   15
##   complement    5
```

When the input is matrix and it contains elements that do not belong to any of
the set, these elements are treated as complement set.


```r
x = list_to_matrix(lt, universal_set = letters)
m = make_comb_mat(x)
m
```

```
## A combination matrix with 3 sets and 8 combinations.
##   ranges of combination set size: c(1, 8).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
##   x      100    1
##     x    010    3
##       x  001    8
##          000    6
## 
## Sets are:
##          set size
##            a    5
##            b   10
##            c   15
##   complement    6
```

The universal set also works for sets as genomic regions.

## Utility functions {#upset-utility-functions}

`make_comb_mat()` returns a matrix, also in `comb_mat` class. There are some
utility functions that can be applied to this `comb_mat` object:

- `set_name()`: The set names.
- `comb_name()`: The combination set names. The names of the combination sets
  are formatted as a string of binary bits. E.g. for three sets of **A**,
  **B**, **C**, the combination set with name "101" corresponds to selecting set
  **A**, not selecting set **B** and selecting set **C**.
- `set_size()`: The set sizes.
- `comb_size()`: The combination set sizes.
- `comb_degree()`: The degree for a combination set is the number of sets that
  are selected.
- `t()`: Transpose the combination matrix. By default `make_comb_mat()`
  generates a matrix where sets are on rows and combination sets are on
  columns, and so are they on the UpSet plots. By transposing the combination
  matrix, the position of sets and combination sets can be swtiched on the
  UpSet plot.
- `extract_comb()`: Extract the elements in a specified combination set. The
  usage will be explained later.

Quick examples are:


```r
m = make_comb_mat(lt)
set_name(m)
```

```
## [1] "a" "b" "c"
```

```r
comb_name(m)
```

```
## [1] "111" "110" "101" "011" "100" "010" "001"
```

```r
set_size(m)
```

```
##  a  b  c 
##  5 10 15
```

```r
comb_size(m)
```

```
## 111 110 101 011 100 010 001 
##   2   1   1   4   1   3   8
```

```r
comb_degree(m)
```

```
## 111 110 101 011 100 010 001 
##   3   2   2   2   1   1   1
```

```r
t(m)
```

```
## A combination matrix with 3 sets and 7 combinations.
##   ranges of combination set size: c(1, 8).
##   mode for the combination size: distinct.
##   sets are on columns
## 
## Combination sets are:
##   a b c code size
##   x x x  111    2
##   x x    110    1
##   x   x  101    1
##     x x  011    4
##   x      100    1
##     x    010    3
##       x  001    8
## 
## Sets are:
##   set size
##     a    5
##     b   10
##     c   15
```

For using `extract_comb()`, the valid combination set name should be from
`comb_name()`. Note the elements in the combination sets depends on the "mode"
set in `make_comb_mat()`.


```r
extract_comb(m, "101")
```

```
## [1] "j"
```

Next we demonstrate a second example, where the sets are genomic regions.
**When the sets are genomic regions, the size is calculated as the sum of the
width of regions in each set (or in other words, the total number of base
pairs).**


```r
library(circlize)
library(GenomicRanges)
lt2 = lapply(1:4, function(i) generateRandomBed())
lt2 = lapply(lt2, function(df) GRanges(seqnames = df[, 1], 
	ranges = IRanges(df[, 2], df[, 3])))
names(lt2) = letters[1:4]
m = make_comb_mat(lt2)
set_size(m)
```

```
##          a          b          c          d 
## 1566783009 1535968265 1560549760 1552480645
```

```r
comb_size(m)
```

```
##      1111      1110      1101      1011      0111      1100      1010      1001 
## 197341532 197137160 194569926 198735008 191312455 192109618 192670258 194462988 
##      0110      0101      0011      1000      0100      0010      0001 
## 191359036 184941701 199900416 199756519 187196837 192093895 191216619
```

And now `extract_comb()` returns genomic regions that are in the corresponding combination set.


```r
extract_comb(m, "1010")
```

```
## GRanges object with 5063 ranges and 0 metadata columns:
##          seqnames            ranges strand
##             <Rle>         <IRanges>  <Rle>
##      [1]     chr1     255644-258083      *
##      [2]     chr1     306114-308971      *
##      [3]     chr1   1267493-1360170      *
##      [4]     chr1   2661311-2665736      *
##      [5]     chr1   3020553-3030645      *
##      ...      ...               ...    ...
##   [5059]     chrY 56286079-56286864      *
##   [5060]     chrY 57049541-57078332      *
##   [5061]     chrY 58691055-58699756      *
##   [5062]     chrY 58705675-58716954      *
##   [5063]     chrY 58765097-58776696      *
##   -------
##   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

With `comb_size()` and `comb_degree()`, we can filter the combination matrix as:


```r
m = make_comb_mat(lt)
# combination set size >= 4
m[comb_size(m) >= 4]
```

```
## A combination matrix with 3 sets and 2 combinations.
##   ranges of combination set size: c(4, 8).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##     x x  011    4
##       x  001    8
## 
## Sets are:
##   set size
##     a    5
##     b   10
##     c   15
```

```r
# combination set degree == 2
m[comb_degree(m) == 2]
```

```
## A combination matrix with 3 sets and 3 combinations.
##   ranges of combination set size: c(1, 4).
##   mode for the combination size: distinct.
##   sets are on rows.
## 
## Combination sets are:
##   a b c code size
##   x x    110    1
##   x   x  101    1
##     x x  011    4
## 
## Sets are:
##   set size
##     a    5
##     b   10
##     c   15
```

For the complement set, the name for this special combination set is only
composed of zeros.


```r
m2 = make_comb_mat(lt, universal_set = letters)
comb_name(m2) # see the first element
```

```
## [1] "111" "110" "101" "011" "100" "010" "001" "000"
```

```r
comb_degree(m2)
```

```
## 111 110 101 011 100 010 001 000 
##   3   2   2   2   1   1   1   0
```

If `universal_set` was set in `make_comb_mat()`, `extract_comb()` can be
applied to the complement set.


```r
m2 = make_comb_mat(lt, universal_set = letters)
extract_comb(m2, "000")
```

```
## [1] "a" "b" "f" "p" "u" "z"
```

```r
m2 = make_comb_mat(lt, universal_set = letters[1:10])
extract_comb(m2, "000")
```

```
## [1] "a" "b" "f"
```

When `universal_set` was set, `extract_comb()` also works for genomic region sets.

## Make the plot {#upset-making-the-plot}

Making the UpSet plot is very straightforward that users just send the
combination matrix to `UpSet()` function:


```r
UpSet(m)
```

<img src="08-upset_files/figure-html/unnamed-chunk-22-1.png" width="480" style="display: block; margin: auto;" />

By default the sets are ordered by the size and the combination sets are
ordered by the degree (number of sets that are selected).

The order is controlled by `set_order` and `comb_order`:


```r
UpSet(m, set_order = c("a", "b", "c"), comb_order = order(comb_size(m)))
```

<img src="08-upset_files/figure-html/unnamed-chunk-23-1.png" width="480" style="display: block; margin: auto;" />

Color of dots, size of dots and line width of the segments are controlled by
`pt_size`, `comb_col` and `lwd`. `comb_col` should be a vector corresponding
to the combination sets. In following code, since `comb_degree(m)` returns a
vector of integers, we just use it as index for the color vector.


```r
UpSet(m, pt_size = unit(5, "mm"), lwd = 3,
	comb_col = c("red", "blue", "black")[comb_degree(m)])
```

<img src="08-upset_files/figure-html/unnamed-chunk-24-1.png" width="480" style="display: block; margin: auto;" />

Colors for the background (the rectangles and the dots representing the set is
not selected) are controlled by `bg_col`, `bg_pt_col`. The length of `bg_col`
can have length of one or two.


```r
UpSet(m, comb_col = "#0000FF", bg_col = "#F0F0FF", bg_pt_col = "#CCCCFF")
```

<img src="08-upset_files/figure-html/unnamed-chunk-25-1.png" width="480" style="display: block; margin: auto;" />

```r
UpSet(m, comb_col = "#0000FF", bg_col = c("#F0F0FF", "#FFF0F0"), bg_pt_col = "#CCCCFF")
```

<img src="08-upset_files/figure-html/unnamed-chunk-25-2.png" width="480" style="display: block; margin: auto;" />

Transposing the combination matrix swtiches the sets to columns and
combination sets to rows.


```r
UpSet(t(m))
```

<img src="08-upset_files/figure-html/unnamed-chunk-26-1.png" width="288" style="display: block; margin: auto;" />

As we have introduced, if do subsetting on the combination sets, the subset of
the matrix can be visualized as well:


```r
UpSet(m[comb_size(m) >= 4])
UpSet(m[comb_degree(m) == 2])
```

<img src="08-upset_files/figure-html/unnamed-chunk-28-1.png" width="672" style="display: block; margin: auto;" />

Following compares the different mode in `make_comb_mat()`:


```r
m1 = make_comb_mat(lt) # the default mode is `distinct`
m2 = make_comb_mat(lt, mode = "intersect")
m3 = make_comb_mat(lt, mode = "union")
UpSet(m1)
UpSet(m2)
UpSet(m3)
```

<img src="08-upset_files/figure-html/unnamed-chunk-30-1.png" width="480" style="display: block; margin: auto;" />

For the plot containing complement set, there is one additional column showing
this complement set does not overlap to any of the sets (all dots are in grey).


```r
m2 = make_comb_mat(lt, universal_set = letters)
UpSet(m2)
```

<img src="08-upset_files/figure-html/unnamed-chunk-31-1.png" width="480" style="display: block; margin: auto;" />

If you already know the size for the complement set, you can directly assign it 
by `complement_size` argument in `make_comb_mat()`.


```r
m2 = make_comb_mat(lt, complement_size = 10)
UpSet(m2)
```

<img src="08-upset_files/figure-html/unnamed-chunk-32-1.png" width="480" style="display: block; margin: auto;" />

For the case where the universal set is smaller than the union of all sets.


```r
m2 = make_comb_mat(lt, universal_set = letters[1:10])
UpSet(m2)
```

<img src="08-upset_files/figure-html/unnamed-chunk-33-1.png" width="480" style="display: block; margin: auto;" />

By default the empty combination sets are removed from the plot, but they can be 
kept by setting `remove_empty_comb_set = FALSE` in `make_comb_mat()`.


































































