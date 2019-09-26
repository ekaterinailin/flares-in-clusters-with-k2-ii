

# Membership matching

## Content

- Overview
- How to work through this section / structure
- Membership studies
    - M67
    - Hyades
    - Praesepe
    - Pleiades
    - M35
    - Ruprecht 147
- Random notes, maybe interesting for further work

## Overview

This document explains how different open cluster membership studies are merged into a final cluster sample. We used the mean membership probability, and required the targets in our final sample to have membership probability $p>0.8$. If a membership study produced binary information on membership only, or if a meta study on membership was conducted in a source we included, we detailed how we assigned membership probabilities here.

The membership studies were presented in heterogeneous data formats. Notebooks starting with a `0x_` contain the preparatory steps to producing merge-able data from the archive material. Notebooks starting with a `1x_` do the bulk of the work. `2x_` do the post-analysis plotting. Intermediate and final results are stored in the folder structure explained below.

## How to work through this section / structure

**If you wish to reproduce the results from this section run the notebooks in order of numbering and follow any instructions given in the cells. You can skip preparatory notebooks to save time in the first iteration. You have to install all necessary packages, and please run the tests before you begin.**

Like this:

    pip install -r requirements.txt

    AND EITHER

    pytest catalog_matching

    OR

    python3 -m pytest catalog_matching

### Structure

1. All preparatory notebooks have a prefix starting with `0x_`. Their results are saved in the `downloads/` folder.

2. All main body cluster matching notebooks have a prefix starting with `1x_`. There is only one such notebook in this section. Just running the entire notebook will produce the final membership tables for all six clusters in this work. For the matching with the K2 archive follow the steps in the notebook one by one. You can also just reproduce results for an individual cluster if you execute the first cells and then the corresponding cluster section.

3. All subsequent analysis and plotting notebooks have a prefix starting with `2x_`. They are used to produce plots that appear in Ilin+(in prep.).

4. Folder `downloads/`
    - archive membership tables as download from Vizier or paper supplements and their reformatted versions, sorted into subfolders by `clusters_<author>/`.
    - `superstamps/` with tabled positions and identifiers for PSF derived light curves for the mre distant clusters M35 and M67
    - a table that described the data formats of archive material that is needed in notebook `01_`

5. Folder `catalog_matching/` contains
    - the functions, scripts, and tests used in the notebooks in this section
    - `requirements.txt` with the most important packages used and their versions (install them all with `$pip install -r requirements.txt`)
    - `matched_catalogs/`
        - `membership_matches/`
            - with all likely members ($p>0.8$) as derived from the different sources as `<cluster>_allmembers.csv`
            - `radec/` with the RA and Dec columns to feed the [K2 Search interface](https://archive.stsci.edu/k2/data_search/search.php) (see further instruction on how to do this in notebook `01_`)

        - `k2_matches_to_union/`
            - the output from the K2 search to use in notebook `01_` in `<cluster>_k2_search.txt`.
            - **the resulting tables to use in the remainder of the project _Flares in Clusters II_ in `<cluster>_k2_mem.csv`.**
            - `members.txt` protocols the number of members with $p>0.8$, how many of these are equipped with short and long cadence light curves from K2. Superstamp derived LCs for M35 and M67 are not included in that table.

6. Folder `plots/` contains
    - `final_sample/`
        - membership histograms for each cluster, broken down by source
        - a pickled dictionary required to produce them.
    - `matched_catalogs/`
        - histograms of mutual distances between matched targets in different studies. The sources are either individual studies or merged ones. In the latter case, this is signified by a composite name in the legend. The filename is read as
            `<cluster>_<study#1>_<study#2>_..._<study#n>_maxradius<cone search radius>arcsec.png`

In the following, we detailed the individual studies included in our work, how we assigned probabilities, other useful notes, and remarks.

## Membership studies

We set the threshold membership probability for a target in our sample to $p=0.8$. Meta studies, and membership studies with classifiers as results were assigned membership probabilities as detailed below.


The following membership studies were used:


study/mission     | source         | membership information quantification   | clusters covered
------------------|----------------|-----------------------------------------|---------------------
Gaia DR2          | Gaia Collaboration 2018    | members list                      | Hyades, M35, Rup 147, Pleiades, Praesepe
DANCe             | Bouy+2015           | probability                             | M35
Gaia DR2          | Cantat-Gaudin+2018  | probability                             | M35, Rup 147, Pleiades, Praesepe
various                 | Curtis+2013         | classifier                             | Rup147
various            | Douglas+2014   | probability                             | Hyades, Praesepe
various               | Douglas+2017   | meta study, members list                | Praesepe
Gaia DR2          | Gao 2018       | probability > .65                       | M67
various             | Gonzalez 2016   | classifier                            | M67
DANCe             | Olivares+2018      | probability                             | Pleiades, Rup 147
DANCe             | Olivares+2019      | probability                             | Rup 147
various             | Rebull+2017         | meta study, classifier                 | Praesepe
various             | Rebull+2016        | meta study, members list                    | Pleiades
Gaia DR1          | Reino+2018          | probability                             | Hyades


Targets that appear in the Gaia HRD diagrams ([Gaia Collaboration 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..10G/abstract), [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/616/A10)) increase our confidence about their membership if it less certain in some other catalog. All open clusters except for M67 appear in that study. If they do not appear in the diagrams we do take this for a reason to be less confident about a member that would be included in the sample judging by membership probabilities in the other catalogs, because the Gaia selection cuts out strong variables and is incomplete for very faint and very bright targets.

**Meta studies and classification studies are assigned membership probabilities as follows:**

Classification in Gonzalez 2016:

    M = 	Member                    .9
    SM = 	Single member            .9
    BM = 	Binary member            .9
    N = 	Non-member                .1
    SN = 	Single non-member        .1
    BN = 	Binary non-member        .1
    U = 	Unknown membership        .5

Classification in Curtis+2013:

    Y =	Yes, highest confidence member;    .9
    P =	Possible/Probable member;          .7
    N =	Not likely/Non-member;             .1
    B =	photometry consistent with Blue stragglers. .0

Rebull+2016, Douglas17, Gaia Collaboration 2018:

    mem : .9

Classification Rebull+2017 (see notebook `03\_`):

    best : .9
    ok : .6
    else : .1

Multiplicity information was available for:

 - M67 (Gonzalez 2016, EB = eclipsing binary)
 - Pleiades (Olivares+2018, equal mass binary)
 - Praesepe (Douglas+17,Douglas+14)
 - Hyades (Douglas+14)
( - Rup 147 (3 individual targets are EBs))


The sections below list individual sources, add remarks, and potentially useful quotes.

#### M67

- [Gao 2018, Gaia ML membership study](https://ui.adsabs.harvard.edu/abs/2018ApJ...869....9G/abstract)
- Gonzalez 2016 was used in Ilin+2019. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/459/1060)

Nardiello+2016 published [PSF de-trended light curves](groups.dfa.unipd.it/ESPG/Kepler-K2.html) along with their [study](https://academic.oup.com/mnras/article-abstract/463/2/1831/2892841?redirectedFrom=fulltext) of variable stars and exoplanet candidates in M67.

#### Hyades

 In the case of Hyades ($d\approx46$ pc), we propagated the sky coordinates from the studies based on Gaia DR1 and DR2 from epoch 2015.5 to epoch 2000. This ensured that a cross matching radius of 3 arcsec would suffice to match all three tables including high proper motion stars (see notebook `04_`). This was not necessary for the other clusters. We plotted the mutual distances obtained from the cross-matching for each combination of sources and stored them in the `plots/` folder.

- Gaia Collaboration 2018
- Reino+2018 (251 members in Gaia DR1). Do not provide RA and Dec, but HIP, Tycho and Gaia IDs. See preparatory notebook `02_`
- Douglas+2014: [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/ApJ/795/161/table5)

ADQL query example for matching Gaia IDs to their RA and Dec:

    SELECT gaiadr1.gaia_source.ra, gaiadr1.gaia_source.dec, gaiadr1.gaia_source.source_id
    FROM user_eilin.hyades_reino
    JOIN gaiadr1.gaia_source
    ON gaiadr1.gaia_source.source_id = user_eilin.hyades_reino.col1


#### Praesepe

- Gaia Collaboration 2018
- [Cantat-Gaudin+2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..93C/abstract). [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93)
- Douglas+2014 [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/ApJ/795/161/table5). Original version of target list included membership probabilities from Klaus and Hillenbrand 2007, and summarized in Agueros+2011 (see below) Membership probabilities given in percent.
- Douglas+2017 use the same sample as Douglas+2014 but with membership probabilties >70%. Douglas+2014 was used in Ilin+2019. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/ApJ/842/83/table3) No individual membership probabilities given.
- Rebull+2017 find a few dozens fewer members than Douglas, [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJ/839/92) No membership probabilities. They retain _likely_ members from other studies, see quotes below.

_Rebull+2017 Section 2.5.2.: Target List_

_[Merged] membership lists in Klein-Wassink (1927), Jones & Cudworth (1983), and Jones & Stauffer (1991), with some candidate members deleted due to discrepant photometry or radial velocities (RVs). This list was then merged with half a dozen recent proper motion membership studies (Adams et al. 2002; Kraus & Hillenbrand 2007; Baker et al. 2010; Boudreault et al. 2012; Khalaj & Baumgardt 2013; Wang et al. 2014), retaining stars considered as likely members in those papers. We then merged this Praesepe membership catalog with the list of all stars observed in K2 Campaign 5 within programs targeting Praesepe. About 600 did not have K2 LCs, sometimes due to the star falling in CCD gaps or just completely outside the K2 FOV; in other cases, the star may have been observable, but no LC was obtained. At this point, then, we have a set of 984 Praesepe members or candidate members with K2 LCs._

_Rebull+2017 Section 2.5.4: Membership analysis_

_From the initial sample of 984 LCs, then, 943 are members of Praesepe by our criteria._

Rebull+2017 Section 2.6: Bright and faint limits

_We have dropped the brightest ( Ks < 6) and retained the rest. There are two stars with Ks < 6, one of which we determined to be periodic, and both of which are listed in Appendix F. Both of them are discarded from our sample as too bright. There are 21 stars with 6 < Ks < 8; K s =8 is roughly an F5 spectral type. At least 11 of them are likely pulsators (with 6 more likely pulsators that have fainter K s ); see Appendix D. We have left these in the sample to allow for comparison to ourPleiades work (which also included likely pulsators), but have identified those pulsators where necessary in the remaining discussion._

_Agueros+2011 Section 2, as cited as PaperI in Douglas 2014_

_Kraus & Hillenbrand (2007) combined data from SDSS (York et al. 2000), the Two Micron All Sky Survey (2MASS; Skrutskie et al. 2006), and USNO-B1.0 (Monet et al. 2003) to calculate proper motions and photometry for several million sources within 7 deg of Praesepe’s center. This census covers a larger area of sky and is deeper than any previous proper motion study of the cluster. The resulting catalog includes 1129 candidate members with membership probability >50% (hereafter referred to as the P50 stars); 442 were identified as high-probability candidates for the first time. Kraus & Hillenbrand (2007) estimated  hat their survey is >90% complete across a wide range of spectral types, from F0 to M5. 16_

#### Pleiades

- Gaia Collaboration 2018
- [Cantat-Gaudin+2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..93C/abstract). [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93)
- Rebull+2016 (table 2 and table 3) was used in Ilin+2019. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/152/113). No membership probabilities given, literature research is detailed in the paper.
- Olivares+2018 DANCe Gaia DR2 survey [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/617/A15)

_Rebull+2016 Section 2.5: Membership and Definition of Sample_

_In order to establish the best possible set of Pleiades members, we evaluated each object using a combination of proper motions and photometric position in an optical color–magnitude diagram (CMD). We primarily used membership probabilities based on recent proper-motion studies (Bouy et al. 2015; see also Sarro et al. 2014 and Lodieu et al. 2012). For objects where the membership probability and the photometric position were inconsistent, we evaluated stars on a case-by-case basis, comparing information from many sources, such as positions and proper motions, radial velocities, X-ray flux, IR flux, and Hα equivalent width. These values are from the literature [...]_

_Olivares+2018: Gaia DR2 Bayesian model_

_From simulations, we estimated that this statistical tool renders kinematic (proper motion) and photometric (luminosity) distributions of the cluster population with a contamination rate of 5.8 ± 0.2%. The luminosity distributions and present day mass function agree with the ones found in a recent study, on the completeness interval of the survey. At the probability threshold of maximum accuracy, the classifier recovers ≈90% of the recently published candidate members and finds 10% of new ones._


#### M35

- Gaia Collaboration 2018 (the authors note that M35 is affected by interstellar reddening)
- [Cantat-Gaudin+2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..93C/abstract). [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93)
- Bouy+2015 DANCe study. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/575/A120)

[Soares-Furtado+2017](https://k2.hatsurveys.org/archive/) extracted PSF derived light curves from K2 Campaign 0 superstamps. This is the corresponding [paper](https://ui.adsabs.harvard.edu/link_gateway/2017PASP..129d4501S/doi:10.1088/1538-3873/aa5c7c).


#### Ruprecht 147

- Gaia Collaboration 2018
- - [Cantat-Gaudin+2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..93C/abstract). [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/618/A93)
- Curtis+2013 present a very detailed membership discussion, focusing on individual members. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/145/134). Providing 2MASS IDs, not coordinates. Membership probability combines:
    - proper motion radial distance from cluster value
    - Lick/Palomar radial velocities
    - Hectochelle radial velocities
    - 2MASS (J-Ks) color-magnitude diagram
    - CFHT/MegaCam
- Olivares+2019 DANCe study. [Source Vizier](http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/625/A115)


Bragaglia+2018 confirm the cluster's old age by the chemical composition of 21 MS and giant members.

## Random notes, maybe interesting for further work

- Hillenbrandt+2018 identified wide binaries in Pleiades and Praesepe.
- A cross-check with _Membership determination of open clusters based on a spectral clustering method_ by Xin-Hua Gao 2018 could be done for some clusters.
- More proper motion studies (esp. Pleiades): Sarro et al. 2014; Lodieu et al. 2012.
- [Yen+2018](https://ui.adsabs.harvard.edu/#abs/2018A&A...615A..12Y/abstract) also determined cluster parameters for Praesepe, Pleiades, Ruprecht 147, and 21 more, including **distance** and **age**.
