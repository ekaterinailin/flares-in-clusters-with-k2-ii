# Stellar Parameters

## Structure

1. All preparatory notebooks have a prefix starting with `0x_`.
    - Before you run any notebooks you can produce the input for `01_` following the recipe in section _TOPCAT cross-matching against 2MASS, PanSTARRS1, and Gaia DR2_ below.
    - `01_ADD_STARHORSE_query_starhorse_for_Teff_by_Gaia_source_id.ipynb` Run this notebook to query GAIA@AIP for StarHorse effective temperatures and calculate average undertainties. Produces tables in `input/StarHorse_output/`.

2. All main notebooks for stellar parameter estimation have a prefix starting with `1x_`.
    - `10_STELLAR_PARAMETERS_Interactive_CMD.ipynb` Run this notebook once for each cluster. You will have to upload the output from your CMD outlier selection to `cmd/outliers/`. More details in section _Interactive CMD selection_. Produces tables in `cmd/cmd_flagged/`
    - `11_STELLAR_PARAMETERS_all_OCs_with_plots.ipynb` Run this notebook once to produce all color-Teff and Teff-R plots, the `luminosities/` and `altai/` tables

3. All subsequent analysis and plotting notebooks have a prefix starting with `2x_`. They are used to produce plots and/or tables that can appear in Ilin+(in prep.), and were not produced in `1x_`.
    - `21_SAMPLE_DESCRIPTION.ipynb`

4. Folder `input/`
    - `members_input/`
        - `<cluster>_k2_mem.csv` the result from membership matching (see 1st part of the project in `Membership_Matching/`
    - `TOPCAT_output/`
        - `<cluster>_topcat_<XX>.csv` results from matching `<cluster>_k2_mem.csv` remotely to Gaia DR2, PanSTARRS1, and 2MASS databases. `<XX>` is either `K2`, `SF`, `N`, or `w_duplicates`.
        - `<cluster>_topcat.csv` refers to merged and cleaned versions of `<cluster>_topcat_<XX>.csv` for Hyades, M35, and M67. Details in notebook `01_` in this folder. Otherwise it is just the TOPCAT cross-matching result.

        _The following tables exist for all clusters except Hyades_

        - `<cluster>_gaiaarchive.csv` was produced in notebook `02_` in this folder, and uploaded to Gaia ADQL interface to obtain `e_bp_min_rp_percentile_lower` and `_upper` that cannot be queried from the TOPCAT interface.
        - `<cluster>_gaiaarchive_results.csv` is the respective output from the cross-matching against the Gaia archive
        - `<cluster>_topcat_w_extinction_percentiles.csv` is `<cluster>_topcat.csv` with two additional columns with uncertainties on extinction

    __The final tables from gathering astrometric and photometric data from Gaia, PanSTARRS, and 2MASS are `hyades_topcat.csv`, and `<cluster>_topcat_w_extinction_percentiles.csv`. These make the basis for the rest of the project part.__

5. Folder `cmd/`
    - `outliers/`
        - `<cluster>_<band1>_<band2>.txt` lists of target IDs, two for each cluster. These target are outliers in the band1-band2 CMDs.
    - `cmd_flagged/`
        - `<cluster>_cmd.csv` same as teh final output from 4., but now without the targets that do have light curves, and with an extra column called `todrop` that indicates if these targets fall off the MS. **These are the tables to proceed with.**

6. Folder `ancillary/` contains the log files with info on when code was run and some additional flags.

6. Folder `clusters/`

    - `Teff_bins_merged.csv` correspondence of SpT and Teff taken from Pecaut and Mamajek (2013)
    - `cluster_parameters_merged.csv` cluster parameters as detailed in Appendix Table B.1
    - `open_clusters_within_300pc.csv` preliminary list of Open Clusters within 300 pc of the Sun

6. Folder `altai/` contains lists of EPIC IDs, and Soares-Furtado et al. (2017) IDs for M35 (NGC 2168) or Nardiello et al. (2016) IDs for M67

7. Folder `luminosities/` contains a table for each cluster with all determined stellar parameters.

8. Folder `flares/` contains the flare tables obtained from AltaiPony

11. Folder `plots/` contains
    - `<cluster>_Teff_spread_all.png` Detailed plot of all Teffs used from different color-temperature relations, StarHorse, and Gaia Apsis
    - `<cluster>_Teff_R.png` Teff vs. stellar radius plot with uncertainties for each cluster

## Introduction

In this part, we started with the membership tables from `Membership_Matching/` for the six clusters.  We then proceeded to cross-match remotely against the 2MASS, PanSTARRS1 (PS1) and Gaia DR2 data bases to to obtain, JHK, grizy, and Gaia photometry. We then flagged obvious non-main sequence stars and outliers in CMDs that we derived for each cluster as outliers (evolved stars, WDs, ...). We transformed grizy from PS1 to SDSS, applied quality cuts, and corrected for extinction using 3D dust maps, or internal corrections in the case of Gaia $B_P-R_P$. We then used color-temperature relations (CTR) to find effective temperatures using empirical relations from multiple studies, inclusing Boyajian et al. (2013) and Mann et al. (2015). From temeperatures we derived stellar radii, again, using empirical relations, this time from Yee et al. (2017) (and references therein), and Mann et al. (2015). Using empirical spectra from Yee et al. (2017) and model spectra from Sarah J. Schmidt we determined the quiescent luminosity of the stars in our sample in the Kepler band.

We also compiled ages and metallicities from various studies conducted in the last 20 years. [Fe/H] values were used, for instance, as parameters in CTRs.

## TOPCAT cross-matching against 2MASS, PanSTARRS1, and Gaia DR2

*Recipe:*

1. upload input table from `input/members_input/`.
2. click: sky cross-match against remote tables [...]
3. type: radius 2. arcsec
4. use RAJ2000_X and DEJ2000_X with X being the suffixes as coordinates
5. Select PanSTARRS DR1, Gaia DR2, or 2MASS
6. use _All_ for column renaming and \_PS1, \_2MASS, or \_Gaia as suffixes, find mode: _Each_\*
7. Go!
8. Save resulting table as .csv file to `input/TOPCAT_output/<clustername>_topcat_<XX>.csv`, use suffixes `<XX>` only if needed.


##### *How many entries did we find for each cluster and database?*

|matches to  | PS1  | 2MASS | Gaia     | upload| suffix\*\*
|------------|------|-------|----------|-------|
|Hyades      | 262  | 279   | 278/538  | 722   | K2\*\*\*
|M35         |5/156 | 5/155 | 5/156  | 1614  | K2/SF
|Ruprecht 147| 96   | 97    | 97  | 238   | K2
|Pleiades    | 938  | 942   | 1240  | 939   | K2
|Praesepe    | 1163 | 1148  | 1157 | 1157  | K2
|M67         | 502/1018 | 496/1006  | 497/1016 | 2020  | K2/N

\*Find mode should be _Each_ (see TOPCAT help):

    Def.: Load a new table with the same number of rows as the local table, in the same order. The best remote match, if any, appears alongside the local row, but if there is no match the remote columns are blank.

\*\* If suffix was only K2, no suffix was added at all.

\*\*\*For the Hyades we matched both the `_GaiaC` and `_K2` coordinates with Gaia DR2.


### Documentation of photometry sources

- [Gaia DR2 main table docs](https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html)
- [2MASS Vizier docs](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/281/2mass6x&-out.add=.)
- [PanSTARRS DR1 Vizier docs](http://vizier.u-strasbg.fr/viz-bin/VizieR-2), [PanSTARRS general page](https://panstarrs.stsci.edu/), [PS1 Glossary](https://outerspace.stsci.edu/display/PANSTARRS/PS1+Glossary+of+PS1-specific+terms)

In the following, we present some details on quality cuts that we applied to PanSTARRS, 2MASS, and Gaia photometry, and quote relevant documentation bits.

#### Gaia DR2

    phot_g_mean_mag : G-band mean magnitude (float, Magnitude[mag])

    Mean magnitude in the G band. This is computed from the G-band mean flux applying the magnitude zero-point in the Vega scale.

    No error is provided for this quantity as the error distribution is only symmetric in flux space. This converts to an asymmetric error distribution in magnitude space which cannot be represented by a single error value.

We computed flux over flux_error for Gaia G, BP, and RP instead. It should be >=10. We used Gaia DR2 internal extinction correction E(BP-RP) using the following columns from the main archive:

    e_bp_min_rp_val : line-of-sight reddening E(BP-RP) (float, Magnitude[mag]). Estimate of redenning E[BP-RP] from Apsis-Priam.

    e_bp_min_rp_percentile_lower / _upper : 16th and 84 percentile of E(BP-RP)

Derivation of these extinction values can be found in [Gaia DR2 documentation](https://gea.esac.esa.int/archive/documentation/GDR2/Data_analysis/chap_cu8par/sec_cu8par_process/ssec_cu8par_process_priamextinction.html)

We used the mean between upper and lower percentile for a simple estimate of the final uncertainty on $BP-RP$.

#### 2MASS Vizier

    Note (1)  : Photometric quality flag (ph_qual). Three character flag, one character per band [JHK], that provides a summary of the net quality of the default photometry in each band, as derived from the Read Flag (rd_flg), measurement uncertainties ([jhk]_msig), scan signal-to-noise ratios ([jhk]_snr), frame-detection statistics (ndet), and profile-fit reduced chi-squared values ([jhk]_psfchi). The value for ph_qual is set for a band according to the precedence of the table below. For example, a source that is tested and meets the conditions for category "X" is not tested for subsequent qualities.
    [...]
    A = 	Detections in any brightness regime where valid measurements were made (rd_flg="1","2" or "3") with [jhk]_snr>10 AND [jhk]_msig<0.10857.
    [...]


&rArr; We require the measurements for each band to be A.

#### PanSTARRS


There are multiple flags explained [here](http://cluster.shao.ac.cn/wiki/index.php/Ps1-flags):

     ==objflags
    Value	Name	Description
    1	QF_OBJ_EXT	extended in our data (eg, PS)
    2	QF_OBJ_EXT_ALT	extended in external data (eg, 2MASS)
    4	QF_OBJ_GOOD	good-quality measurement in our data (eg,PS)
    8	QF_OBJ_GOOD_ALT	good-quality measurement in external data (eg, 2MASS)
    16	QF_OBJ_GOOD_STACK	good-quality object in the stack (> 1 good stack)
    32	QF_OBJ_SUSPECT_STACK	suspect object in the stack (> 1 good or suspect stack, less tham 2 good)
    64	QF_OBJ_BAD_STACK	good-quality object in the stack (> 1 good stack)

&rArr; We required `QF_OBJ_GOOD` to be set.

##### Conversion from PanSTARRS1 to SDSS photometry.

[Finkbeiner et al. (2016)](http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/meta): SDSS and PS1 band magnitudes are similar but not identical. A corresponding function is implemented in `utils.py` from Table 2 in Finkbeiner et al. (2016)


## Interactive CMD selection

Notebook `10_` contains the interactive CMD generator, and cuts the table of targets to those that have light curves. You'll find the following recipe implemented:

1. Choose the `C` variable to be one of the clusters. Drop all targets without LCs, do the extinction correction, transformation of photometric systems, etc., i.e. all that is implemented in `utils.prepare_and_validate_stars`.
2. Generate CMDs for all colors used in CTRs later\*. Open the `.hmtl` file in a browser (Firefox and Chrome both may work).
3. Select data points that fall off the MS most obviously.
5. Click the _Save_ button and save to `cmd/outliers/<cluster>_<band1>_<band2>.txt`.
6. Go back to the notebook and read in all outlier files, add the corresponding flags to the `todrop` column.
7. Write out the resulting table to `cmd/cmd_flagged/<cluster>_cmd.csv`

Repeat for all clusters.

\* We used g-i, g-J, g-H, g-K, r-z, r-J, BP-RP, g-z. For the CTRs, we did not use BP-RP because results from its CTR showed systematic offset in comparison with spectroscopic studies (Rocio Kiman, priv. comm.). We did not use r-z because it is metallicity dependent, and Boyajian et al. (2013) did not include [Fe/H] as a parameter in their CTRs (Sarah J. Schmidt, priv. comm.).

## Dust extinction

We corrected the photometric data for reddening in JHK and grizy bands using 3D dust maps from [Green et al. 2019](https://arxiv.org/abs/1905.02734). Instruction for usage available from their [webpage](http://argonaut.skymaps.info/). We used the "converged" (`bayestar2019_flag == 3`)and "reliable_dist" (`bayestar2019_flag == 4`) flags to indicate if an extinction value was numerically reliable and lied withing the distance in which extinction could be determined.

We do not apply extinction correction to the Hyades, Pleiades, and Praesepe (see notebook `11_`, `extinction` flag passed to constructor), as they are close enough.

For $B_P-R_P$ we used internal correction from Gaia DR2 (see above).

## Color-Temperature Relations and effective temperatures

In general, empirical relations suffer from systematic errors that stem both from the different methods applied, and from sample selection biases. We therefore used multiple empirical relations in their appropriate ranges to obtain multiple $T_\mathrm{eff}$ estimates from which we then obtained a weighted mean and uncertainty (it is called `Teff_median` and `Teff_std` although this is technically the wrong term).

`cmd.py` contains the relevant functions.


#### Boyajian et al. (2013)

We used Table 8 from [Boyajian et al. (2013)](http://iopscience.iop.org/article/10.1088/0004-637X/771/1/40/meta#apj474530fd2). The CTRs are given as cubic fit in Eqn. (2) in that paper. The authors determine CTRs from a range of interferometrically characterized stars using g-z, g-i, g-r, g-J, g-H, and g-K color from SDSS and Johnson magnitudes for A to K stars.

They use Johnson JHK, not 2MASS JHK, so we did the same transformation as the authors from [2MASS to Bessel-Brett](http://www.astro.caltech.edu/~jmc/2mass/v3/transformations/) ([Carpenter (2001)](https://ui.adsabs.harvard.edu/link_gateway/2001AJ....121.2851C/doi:10.1086/320383)) and from Bessel-Brett to Johnson with [Bessel and Brett (1988)](http://articles.adsabs.harvard.edu/full/1988PASP..100.1134B). This gives us J-H, H-K, and J-K in Johnson system. We can then convert to J, H, and K in the Johnson system directly because J_BB = J_Johnson, because they are both Vega calibrated systems.

We use the CTRs given in Table 8 in the appropriate magnitude ranges, and, if available, only those that are derived from stars with $T_\mathrm{eff} < 6750$ K and that were not only observed in 2MASS (footnotes _c_ and _d_ in Table 8).

#### Mann et al. (2015) (and Erratum)

We used Table 2 found in the Erratum to [Mann et al. (2015)](https://ui.adsabs.harvard.edu/#abs/2015ApJ...804...64M/abstract). The CTRs are given as quadratic polynomial fits to r-z, r-J, or BP-RP with extra information from metallicity or J-H if available. The relation are only valid, if metallicity is sufficiently close to solar, which is satisfied for all of our clusters:

    These relations are valid for 0.1- 0.7 RSun (spectral types K7–M7, 4.6 < MKs < 9.8 , 2700 < Teff < 4100 K) and -0.6 < [Fe/H] < 0.5.

We apply the temperature cuts a posteriori, dropping Teff derived with these CTRs if they are to high or too low in the range given in the quote above, we also apply the cuts on $M_{\mathrm{K}_\mathrm{s}}$.

[Fe/H] is not sampled evenly: later than M4 are mostly metal poor, late K dwarfs are rather metal rich.


##### Dwarfs with SpT later than M5.5

For these UCDs, we use Mann+2015 CTR, although their sample is biased towards metal-poor stars there. As Specmatch-Emp does not cover $T_{eff}<3000$K, we use their $T_{eff}$-radius relation from Table 1, eqn. 3 and 4. Checked the [Erratum](http://adsabs.harvard.edu/abs/2016ApJ...819...87M) here, but there is no difference for eqn. 3 and 4.

$$X = T_{eff}/3500$$

$$R_* = 10.5440 - 33.7546 X + 35.1909 X^2 -11.59280 X^3$$

$$R_* = 16.7700 - 54.3210 X + 57.6627 X^2 −19.69940 X^3 + (1.+0.45650 [Fe/H])$$

#### Gaia DR 2 Apsis

Gaia DR2 published effective temperatures for over 160 million sources ([Gaia Collaboration 2018](https://doi.org/10.1051/0004-6361/201833051)). The typical uncertainty is quoted as 324 K, but it is considerably lower for stars above \~4000 K and below \~6700K, so that we adopt 175 K which is above the quoted root-median-squared error in this Teff range ([Andrae et al. 2018](https://doi.org/10.1051/0004-6361/201732516)).

#### Gaia DR 2 with StarHorse 

## Radii

We determined stellar radii using Mann et al. (2015) Erratum Table 1 relation between $T_\mathrm{eff}$ and $R_\*$ for $2700\,\mathrm{K}<T_\mathrm{eff}<3800\,\mathrm{K}$ (Eqn. (3) and (4)). For stars for which we computed radii from temperatures, we also calculated radii directly from $M_{\mathrm{K}_\mathrm{s}}$ as consistency check (Eqn.(1) and (2)). For all other stars we followed [Yee et al. (2017)](https://iopscience.iop.org/article/10.3847/1538-4357/836/1/77). They compiled a library of low-mass stars with well-determined parameters. The library can be accessed via the Python package [SpecMatch-Emp](https://specmatch-emp.readthedocs.io/en/latest/).

## Luminosities

We determined:

1. `Lum_SED`  Planck curve projected luminosity with SED convolved (SED either from SpecMatch-Emp for $T_\mathrm{eff} > 3500$ K or from Sarah J. Schmidt's model spectra for UCDs available from this [Github repository](https://github.com/sjschmidt/FlareSpecPhot/tree/master/qSED)) and stored here under `FlareAnalysisPipeline/opencluster/static/ML_SEDs/`.
2. `Lum_Kepler`  Kepler projected luminosity with SED: as in 1. but convolved with Kepler response
3. `Lum_SB`  Stefan-Boltzmann luminosity from $T_\mathrm{eff}$ and $R_*$.

## Ages, metallicities, and distances

We compiled a non-comprehensive list of age, distance and metallicty studies for the six clusters in the paper's appendix in Table B.1. There, we also give the values we adopted. For distance, we used Gaia distance when available and good quality, and took the median value from these stars to fill in those without measured Gaia parallaxes.
