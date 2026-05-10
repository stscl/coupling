# coupling

***Coupling** Coordination Analysis*

*coupling* is an R package for coupling coordination analysis based on
coupling coordination degree (CCD) models. It incorporates metacoupling
frameworks to characterize cross-scale linkages, flows, and feedbacks
within and among coupled systems.

> *Refer to the package documentation
> <https://stscl.github.io/coupling/> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=coupling) with:

``` r

install.packages("coupling", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/coupling) with:

``` r

install.packages("coupling",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on
  [GitHub](https://github.com/stscl/coupling) with:

``` r

if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/coupling",
                         build_vignettes = TRUE,
                         dep = TRUE)
```

## References

Tang, P., Huang, J., Zhou, H., Fang, C., Zhan, Y., Huang, W., 2021.
Local and telecoupling coordination degree model of urbanization and the
eco-environment based on RS and GIS: A case study in the Wuhan urban
agglomeration. Sustainable Cities and Society 75, 103405.
<https://doi.org/10.1016/j.scs.2021.103405>.

Wang, S., Kong, W., Ren, L. and ZHI, D., 2021. Research on misuses and
modification of coupling coordination degree model in China. Journal of
Natural Resources, 36, 793-810.
<https://doi.org/10.31497/zrzyxb.20210319>.

Fan, D., Ke, H. and Cao, R., 2024. Modification and improvement of
coupling coordination degree model. Stat. Decis, 40, 41-46.
<https://doi.org/10.13546/j.cnki.tjyjc.2024.22.007>.

Xie, A., Zhang, F., Ding, Y., Chen, J., Yang, P., Peng, G., 2025.
Exploring the Dynamic Local and Tele‐Coupling Coordination Mechanism of
the Ecosystem Services Supply–Demand and Its Driving Forces: Taking
China’s Yangtze River Economic Belt as an Example. Land Degradation &
Development 36, 3178–3193. <https://doi.org/10.1002/ldr.5559>.

Li, Y., Jia, N., Zheng, L., Yin, C., Chen, K., Sun, N., Jiang, A., Wang,
M., Chen, R., Zhou, Z., 2026. A meta-coupling analysis between
three-dimensional urbanization and ecosystem services in China’s urban
agglomerations. Communications Earth & Environment 7.
<https://doi.org/10.1038/s43247-025-03047-w>.

 
