

# Introduction {#introduction}

Complex heatmaps are efficient to visualize associations between different
sources of data sets and reveal potential patterns. Here the
**ComplexHeatmap** package provides a highly flexible way to arrange multiple
heatmaps and supports self-defined annotation graphics.

## General design {#general-design}



A single heatmap is composed of the heatmap body and the heatmap components.
The heatmap body can be split by rows and columns. The heatmap components are
titles, dendrograms, matrix names and heatmap annotations, which are put on
the four sides of the heamap body. The heatmap components are reordered or
split according to the heatmap body.

<img src="01-introduction_files/figure-html/unnamed-chunk-3-1.png" width="672" style="display: block; margin: auto;" />

A heatmap list is concatenation of a list of heatmaps and heatmap annotations.
Surrounding the heatmap list, there are global-level titles and legends.

One important things for the heatmap list is that rows for all heatmaps and
annotations (it is row annotation if the heatmap list is horizontal.) are all
adusted so that the same row in all heatmaps and annotations corresponds to a
same feature.

<img src="01-introduction_files/figure-html/unnamed-chunk-4-1.png" width="384" style="display: block; margin: auto;" />

The heatmaps and annotations (now it is column annotation) can also be
arranged vertically.

<img src="01-introduction_files/figure-html/unnamed-chunk-5-1.png" width="384" style="display: block; margin: auto;" />

And the heatmap list can also be split by rows and by columns.

The **ComplexHeatmap** package is implemented in an object-oriented way. To
describe a heatmap list, there are following classes:

- `Heatmap` class: a single heatmap containing heatmap body, row/column names,
  titles, dendrograms and row/column annotations.
- `HeatmapList` class: a list of heatmaps and heatmap annotations.
- `HeatmapAnnotation` class: defines a list of row annotations and column
  annotations. The heatmap annotations can be components of heatmap, also they
  can be independent as heatmaps.

There are also several internal classes:

- `SingleAnnotation` class: defines a single row annotation or column
  annotation. The `HeatmapAnnotation` object contains a list of
  `SingleAnnotation` objects.
- `ColorMapping` class: mapping from values to colors. The color mappings of
  the main matrix and the annotations are controlled by `ColorMapping` class.
- `AnnotationFunction` class: constructs user-defined annotations. This is the
  base of creating user-defined annotation graphics.

**ComplexHeatmap** is implemented under **grid** system, so users need to know
basic **grid** functionality to get full use of the package.


## A brief description of following chapters {#a-brief-description-of-following-chapters}

- [**A Single Heatmap**](a-single-heatmap.html)

This chapter describes the configurations of a single heatmap.

- [**Heatmap Annotations**](#heatmap-annotations.html)

This chapter describes the concept of the heatmap annotation and demonstrates
how to make simple annotations as well as complex annotations. Also, the
chapter explains the difference between column annotations and row
annotations.

- [**A List of Heatmaps**](a-list-of-heatmaps.html)

This chapter describes how to concatenate a list of heatmaps and annotations
and how adjustment is applied to keep the correspondence of the heatmaps.

- [**Legends**](legends.html)

This chapter describes how to configurate the heatmap legends and annotation
legends, also how to create self-defined legends.

- [**Heatmap Decoration**](heatmap-decoration.html)

This chapter describes methods to add more self-defined graphics to the
heatmaps after the heatmaps are generated.

- [**OncoPrint**](oncoprint.html)

This chapter describes how to make oncoPrints and how to integrate other
functionalities from **ComplexHeatmap** to oncoPrints.

- [**UpSet plot**](upset-plot.html)

This chapter describes how to make enhanced UpSet plots.

- [**Other High-level Plots**](other-high-level-plots.html)

This chapter describes functions implemented in **ComplexHeatmap** for
specific use, e.g. visualizing distributions.

- [**More Examples**](more-examples.html)

More simulated and real-world examples are demonstrated in this chapter.
