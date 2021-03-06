<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicrobiomeAnalysis: A comprehensive R package for managing and analyzing microbiome data

This package has many unique features:

## The Overview of **MicrobiomeAnalysis**

<img src="./man/figures/Schematic.png" title="The Overview of MicrobiomeAnalysis" alt="The Overview of MicrobiomeAnalysis" style="display: block; margin: auto;" />

## :writing_hand: Authors

[Hua Zou](mailto:zouhua1@outlook.com)

## :arrow_double_down: Installation

### Installation of dependencies

Requiring the R-base version more than **4.1.2**. Some of the
dependencies from CRAN or bioconductor, but others, are only obtained
from github. the version of dependencies are as follows:

-   dplyr (>= 1.0.8)
-   phyloseq (>= 1.38.0)
-   magrittr (>= 2.0.2)
-   purrr (>= 0.3.4)
-   MASS (>= 7.3-55)
-   ggplot2 (>= 3.3.5)
-   tibble (>= 3.1.6)
-   rlang (>= 1.0.2)
-   coin (>= 1.4-2)
-   ggtree (>= 3.2.1)
-   tidytree (>= 0.3.9)
-   IRanges (>= 2.28.0)
-   tidyr (>= 1.2.0)
-   patchwork (>= 1.1.1)
-   ggsignif (>= 0.6.3)
-   metagenomeSeq (>= 1.36.0)
-   DESeq2 (>= 1.34.0)
-   edgeR (>= 3.36.0)
-   vegan (>= 2.5-7)
-   limma (>= 3.50.1)
-   caret (>= 6.0-92)
-   scales (>= 1.1.1)
-   yaml (>= 2.3.5)
-   RColorBrewer (>= 1.1-2)
-   ALDEx2 (>= 1.26.0)
-   ape (>= 5.6-2)
-   ade4 (>= 1.7-18)
-   Rtsne (>= 0.15)
-   ANCOMBC (>= 1.4.0)
-   coin (>= 1.4-2)
-   glmnet (>= 4.1-3)
-   BiocGenerics (>= 0.40.0)
-   biomformat (>= 1.22.0)
-   S4Vectors (>= 0.32.3)
-   Biobase (>= 2.54.0)

There are two methods to install the aforementioned packages.

#### dependencies in CRAN & Bioconductor

``` r
install.packages("pacman")

library(pacman)

pacman::p_load(phyloseq, dplyr, tibble, Biobase, ggplot2)
```

#### dependencies in github

``` r
# Step 1: Install devtools
if (!requireNamespace(c("remotes", "devtools"), quietly=TRUE)) {
  install.packages(c("devtools", "remotes"))
}
library(devtools)
#library(remotes)

# Step 2: install corncob package
devtools::install_github("bryandmartin/corncob")
# remotes::install_github("bryandmartin/corncob")
```

#### dependencies from website

Manual download of RAIDA_1.0.tar.gz from
[here](http://cals.arizona.edu/~anling/software/RAIDA_1.0.tar.gz) and
install locally

``` shell
R CMD INSTALL RAIDA_*.tar.gz
```

``` r
install.packages("RAIDA_1.0.tar.gz", repos = NULL, type = "source")
```

Manual download of mbzinb_0.2.tar.gz from
[here](https://github.com/jchen1981/MicrobiomeDDA/blob/master/mbzinb_0.2.tar.gz)
and install locally

``` shell
R CMD INSTALL mbzinb_*.tar.gz
```

``` r
install.packages("mbzinb_0.2.tar.gz", repos = NULL, type = "source")
```

### Installing MicrobiomeAnalysis

Get the released **MicrobiomeAnalysis** packages from
[here](https://github.com/HuaZou/MicrobiomeAnalysis/releases).

``` r
install.packages("MicrobiomeAnalysis*.tar.gz", repos = NULL, type = "source")
```

the development version from github:

``` r
remotes::install_github("HuaZou/MicrobiomeAnalysis")
```

# :book: Vignette

Using the following command and Choosing the `html` for more details.

``` r
utils::browseVignettes(package="MicrobiomeAnalysis")
```

**[the tutorial of Microbiota Data Analysis by
MicrobiomeAnalysis](https://zouhua.top/MicrobiomeAnalysis_book/) made by
bookdown**.

## :sparkling_heart: Contributing

Welcome any contributions or comments, and you can file them
[here](https://github.com/HuaZou/MicrobiomeAnalysis/issues).

## :trophy: Acknowledgement

Thanks all the developers of the methods integrated into
**MicrobiomeAnalysis**. The reference package is
[microbiomeMarker](https://github.com/yiluheihei/microbiomeMarker). The
R codes of differentail analysis methods are almost from
[microbiomeMarker](https://github.com/yiluheihei/microbiomeMarker).

## :eight_pointed_black_star: Citation

Kindly cite by using `citation("MicrobiomeAnalysis")` if you think
**MicrobiomeAnalysis** helps you. Alternative way is Zou H (2022).
*MicrobiomeAnalysis: An R package for analysis and visualization in
metagenomics*. R package version 1.0.1,
\<URL:<https://github.com/HuaZou/MicrobiomeAnalysis/>\>.
