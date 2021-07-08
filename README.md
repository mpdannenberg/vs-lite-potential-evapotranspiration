# vs-lite-potential-evapotranspiration
Effects of potential evapotranspiration (PET) methods on forward models of tree-ring width

## Citation
Dannenberg, M. P. (2021), Modeling tree radial growth in a warming climate: Where, when, and how much do potential evapotranspiration models matter?, Environmental Research Letters, Accepted.

## Abstract
Process-based models of tree-ring width are used both for reconstructing past climates and for projecting changes in growth due to climate change. Since soil moisture observations are unavailable at appropriate spatial and temporal scales, these models generally rely on simple water budgets driven in part by temperature-based potential evapotranspiration (PET) estimates, but the choice of PET model could have large effects on simulated soil moisture, moisture stress, and radial growth. Here, I use four different PET models to drive the VS-Lite model and evaluate the extent to which they differ in both their ability to replicate observed growth variability and their simulated responses to projected 21st century warming. Across more than 1,200 tree-ring width chronologies in the conterminous United States, there were no significant differences among the four PET models in their ability to replicate observed radial growth, but the models differed in their responses to 21st century warming. The temperature-driven empirical PET models (Thornthwaite and Hargreaves) simulated much larger warming-induced increases in PET and decreases in soil moisture than the more physically realistic PET models (Priestley-Taylor and Penman-Monteith). In cooler and more mesic regions with relatively minimal moisture constraints to growth, the models simulated similarly small reductions in growth with increased warming. However, in dry regions, the Thornthwaite- and Hargreaves-driven VS-Lite models simulated an increase in moisture stress roughly double that of the Priestley-Taylor and Penman-Monteith models, which also translated to larger simulated declines in radial growth under warming. While the lack of difference in the models’ ability to replicate observed radial growth variability is an encouraging sign for some applications (e.g., attributing changes in growth to specific climatic drivers), the large differences in model responses to warming suggest that caution is needed when applying the temperature-driven PET models to climatic conditions with large trends in temperature.

## VS-Lite functions
###VS-Lite's functions were modified from version 2.3 of those provided by Susan Tolwinski-Ward (Tolwinski-Ward et al., 2011, 2013, 2014) at https://www.ncdc.noaa.gov/paleo-search/study/9894. These functions include:
Compute_gE.m: Computes the solar energy scalar (unmodified from version 2.3)
estimate_vslite_params_v3.m: Bayesian calibration of VS-Lite parameters (modified from version 2.3 to allow optional switches for other PET methods)
leakybucket_monthly.m: Calculate monthly "leaky bucket" soil moisture without daily substepping (modified from version 2.3 to allow optional switches for other PET methods)
leakybucket_submonthly.m: Calculate monthly "leaky bucket" soil moisture with daily substepping (modified from version 2.3 to allow optional switches for other PET methods)
VSLite_v3.m: Run VS-Lite model (modified from version 2.3 to allow optional switches for other PET methods)

## PET functions
fao_pm.m: Computes PET with the simplified FAO-56 Penman-Monteith method (Allen et al., 1998)
hargreaves.m: Computes PET with the Hargreaves method (Hargreaves & Samani, 1985; Hargreaves & Allen, 2003)
priestley_taylor.m: Computes PET with the Priestley-Taylor method (Priestley & Taylor, 1972)

## Analysis scripts
All scripts, with descriptions, are listed in main.m, which runs them in the correct order.

## References
Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300(9), D05109.

Hargreaves, G. H., & Allen, R. G. (2003). History and evaluation of Hargreaves evapotranspiration equation. Journal of Irrigation and Drainage Engineering, 129(1), 53–63. https://doi.org/10.1061/(ASCE)0733-9437(2003)129:1(53)

Hargreaves, G. H., & Samani, Z. A. (1985). Reference crop evapotranspiration from temperature. Applied Engineering in Agriculture, 1(2), 96–99.

Priestley, C., & Taylor, R. (1972). On the assessment of surface heat flux and evaporation using large scale parameters. Monthly Weather Review, 100, 81–92.

Tolwinski-Ward, S. E., Anchukaitis, K. J., & Evans, M. N. (2013). Bayesian parameter estimation and interpretation for an intermediate model of tree-ring width. Climate of the Past, 9(4), 1481–1493. https://doi.org/10.5194/cp-9-1481-2013

Tolwinski-Ward, S. E., Evans, M. N., Hughes, M. K., & Anchukaitis, K. J. (2011). An efficient forward model of the climate controls on interannual variation in tree-ring width. Climate Dynamics, 36(11–12), 2419–2439. https://doi.org/10.1007/s00382-010-0945-5

Tolwinski-Ward, S. E., Tingley, M. P., Evans, M. N., Hughes, M. K., & Nychka, D. W. (2014). Probabilistic reconstructions of local temperature and soil moisture from tree-ring data with potentially time-varying climatic response. Climate Dynamics, 44(3–4), 791–806. https://doi.org/10.1007/s00382-014-2139-z

