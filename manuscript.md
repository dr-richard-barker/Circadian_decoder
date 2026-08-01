# Spaceflight Alters Plant Circadian Clock Phase: A ChronoGauge Analysis of NASA GeneLab Arabidopsis Transcriptomics

## Authors

Richard Barker¹,*

¹ Department of Botany, University of Wisconsin-Madison, Madison, WI, USA

*Corresponding author: dr.richard-barker@wisc.edu (or dr.richard.barker@gmail.com)

## Abstract

Spaceflight imposes multiple stressors on plants, including microgravity, altered light regimes, and atmospheric changes. The circadian clock regulates approximately one-third of the plant transcriptome, yet whether spaceflight disrupts circadian timing in plants remains poorly understood. We applied ChronoGauge, a deep learning circadian time predictor, to 23 Arabidopsis thaliana spaceflight transcriptomics datasets from NASA GeneLab, encompassing 603 samples across 13 ecotypes, 7 tissue types, and 12 hardware configurations. ChronoGauge predicted circadian time (CT) for each sample from gene expression alone (mean absolute error = 44.0 minutes on held-out validation data), enabling per-study flight-versus-ground phase comparisons using circular statistics. Random-effects meta-analysis across 18 studies with paired flight and ground controls revealed a small but statistically significant phase advance under spaceflight (pooled shift = −0.09 h, 95% CI [−0.18, −0.001], p = 0.046, I² = 86.5%). The effect was driven primarily by root tissue studies (pooled shift = −0.17 h, p < 0.001, I² = 8.9%) and dark-grown samples (pooled shift = −0.20 h, p = 0.024, I² = 91.8%). Four individual studies showed significant phase shifts after Watson-Williams testing, with shifts ranging from −0.17 h to −0.52 h. A forked analysis comparing the 14-study Barker 2023 cohort (Fork A, pooled = −0.09 h, p = 0.019) with all 23 Arabidopsis spaceflight studies (Fork B, pooled = −0.09 h, p = 0.047) demonstrated consistent conclusions regardless of dataset scope. Flight samples showed slightly lower circular variance than ground controls (0.488 vs 0.507), suggesting marginally more coherent clock output under spaceflight. Circadian trajectory analysis using t-SNE on the 100-dimensional sub-model prediction vectors revealed that the multivariate circadian fingerprint correlated with phase shift magnitude (Spearman ρ = 0.633, p = 0.005). Gene set enrichment analysis of the four significant phase-shift studies linked the phase advance in OSD-193 to significant downregulation of the circadian rhythm pathway (GO:0007623, NES = −2.25, FDR = 0.004), with core clock genes (LHY, TOC1, PRR5, CCA1) among the leading edge. This study demonstrates the utility of deep learning-based circadian time prediction for retrospective analysis of cross-sectional transcriptomics data and provides the first systematic quantification of circadian phase disruption in spaceflight-grown plants.

**Keywords:** circadian clock, spaceflight, Arabidopsis thaliana, deep learning, ChronoGauge, NASA GeneLab, transcriptomics, meta-analysis

---

## Introduction

The circadian clock is an endogenous timekeeping mechanism that orchestrates physiological and molecular rhythms with approximately 24-hour periodicity. In plants, the circadian clock regulates roughly one-third of the transcriptome [1], controlling processes including photosynthesis, growth, flowering time, and stress responses [2]. The Arabidopsis circadian oscillator comprises a network of transcriptional-translational feedback loops: the morning complex (CCA1/LHY) represses evening genes (TOC1, GI, ELF4, ELF3, LUX), while sequential PRR proteins (PRR9, PRR7, PRR5) provide daytime repression of CCA1/LHY [3]. Circadian regulation is critical for plant fitness; clock disruption reduces biomass, alters stress responses, and impairs reproductive timing [4].

Spaceflight imposes unique environmental conditions on plants, including microgravity, altered atmospheric composition, modified light quality and quantity, and increased radiation exposure [5]. These conditions may disrupt circadian timing through multiple mechanisms: altered light cues for clock entrainment, mechanical unloading affecting calcium signaling, and changes in metabolic flux. Several studies have reported spaceflight-induced changes in clock gene expression [6,7], and a recent meta-analysis of Arabidopsis spaceflight transcriptomes identified circadian-related gene sets among the most consistently altered pathways [8]. However, these analyses examined differential expression of individual clock genes rather than quantifying circadian phase itself, and a systematic measurement of circadian phase disruption has been lacking.

A major challenge in analyzing circadian rhythms from spaceflight transcriptomics data is that most NASA GeneLab datasets are cross-sectional—samples are harvested at a single time point without circadian time (CT) ground truth. Traditional circadian analysis methods such as JTK_CYCLE [9] and RAIN [10] require time-series data and cannot be applied to single-timepoint samples. Recently, deep learning approaches have emerged that can predict circadian time from a single transcriptomic snapshot. ChronoGauge [11] is a bagging ensemble of 100 neural network sub-predictors trained on Arabidopsis circadian time-course data under constant light, achieving mean absolute errors of approximately 50 minutes on held-out data. Each sub-model uses a distinct feature gene set selected via sequential feature selection, and predictions are aggregated via circular mean. The circular variance across sub-models provides a measure of prediction confidence and clock robustness.

Here, we apply ChronoGauge to the complete collection of Arabidopsis spaceflight transcriptomics data in NASA GeneLab. We perform a forked analysis: Fork A examines the 14 studies from Barker et al. (2023) [8], while Fork B includes all 23 Arabidopsis spaceflight transcriptomics studies with processed expression data. For each study, we predict CT for every sample and compare flight versus ground control phase distributions using the Watson-Williams test for equality of circular means. We then perform random-effects meta-analysis across studies, stratified by tissue, ecotype, hardware, and light regime. This approach enables the first systematic quantification of circadian phase disruption in spaceflight-grown plants.

## Methods

### Data acquisition

We queried the NASA Open Science Data Repository (OSDR) API (https://osdr.nasa.gov) for all Arabidopsis thaliana transcriptomics studies. Fork A comprised 14 studies (GLDS-7, 17, 37, 38, 44, 46, 120, 121, 136, 147, 205, 208, 218, 251; GLDS-213 was unavailable in OSD). Fork B comprised 23 Arabidopsis spaceflight transcriptomics studies. Combined, 26 unique studies were identified. For each study, we downloaded processed expression matrices (DESeq2 median-of-ratios normalized counts for RNA-seq; probeset-level normalized expression for microarrays) and ISA-Tab metadata via the OSDR files API. Three studies were excluded from ChronoGauge analysis: OSD-213 (legacy study not migrated to OSD), OSD-480 (raw FASTQ only, no processed data), and OSD-522 (unnormalized counts only). This yielded 23 studies with usable expression data.

### Metadata curation

ISA-Tab sample files (s_*.txt) were parsed to extract sample-level metadata: condition (spaceflight/ground control), ecotype, tissue, hardware, light regime, and harvest age. Ecotype names were standardized (e.g., "Wasselewskija" mapped to "Ws"). Light regimes were categorized as continuous light, dark, or photoperiod (LD 16:8, LD 12:12, LD 4:2). Studies without a Spaceflight factor (OSD-46, radiation; OSD-136, hypobaria) were identified as ground-only studies and excluded from flight-versus-ground comparisons. OSD-208 was identified as a ground-only root apex study. OSD-251 and OSD-346 had flight samples only (no ground controls) and were excluded from phase-shift comparisons. The final harmonized metadata table comprised 1,169 samples across 26 studies, with 792 flight and 324 ground control samples.

### ChronoGauge application

ChronoGauge [11] is a bagging ensemble of 100 neural network sub-predictors, each trained on a distinct feature gene set (11–32 genes per sub-model) selected via sequential feature selection from Arabidopsis circadian time-course data under constant light. We used pre-trained models from HuggingFace, with separate ensembles for RNA-seq (100 models), ATH1 microarray (100 models), and AraGene microarray (100 models). For each study, we: (1) loaded the normalized expression matrix, (2) mapped gene identifiers to AGI codes (AT*G* format), (3) fit a StandardScaler on ChronoGauge's training data and applied it to the study expression matrix, (4) ran all 100 sub-models with their respective feature gene sets, and (5) aggregated predictions via circular mean. The circular variance across sub-models (1 − R, where R is the mean resultant length) served as a clock robustness metric, with lower values indicating more coherent ensemble predictions. Sub-models were skipped if fewer than 80% of their feature genes were present in the study expression matrix.

### Statistical analysis

For each study with both flight and ground control samples (minimum 3 replicates each), we compared predicted CT distributions using the Watson-Williams test for equality of circular means [12]. Phase shift was computed as the circular difference between flight and ground mean CTs, expressed in the range [−12, +12] hours. Negative values indicate phase advance (flight samples predicted at earlier CT than ground controls). Bootstrap 95% confidence intervals (1,000 iterations) were computed for each phase shift. Random-effects meta-analysis was performed using the REML estimator via the metafor R package [13]. Heterogeneity was assessed with Cochran's Q test and the I² statistic. Stratified meta-analyses were performed by tissue type, hardware configuration, and light regime. All analyses were conducted for both Fork A and Fork B independently.

### Supplementary analyses

Core clock gene expression (CCA1, LHY, TOC1, PRR9, PRR7, PRR5, GI, ELF4, ELF3, LUX, TIC) was extracted from each study's expression matrix, z-scored per gene, and visualized as a heatmap. Circadian fingerprint analysis compared circular variance (clock robustness) between flight and ground conditions. Principal component analysis was performed on the 100-dimensional sub-model prediction vectors to visualize sample clustering by condition and study.

#### Circadian trajectory analysis

To explore whether the multivariate circadian fingerprint captured by ChronoGauge's 100 sub-model predictions reflects the magnitude of phase disruption, we performed dimensionality reduction on the 100-dimensional prediction vectors (603 samples × 100 features). We applied t-distributed stochastic neighbor embedding (t-SNE; perplexity = 30, 1,000 iterations, random seed = 42) and principal component analysis (PCA) after standardizing features with StandardScaler. For each study with paired flight and ground controls, we computed the Euclidean distance between flight and ground centroids in both the t-SNE and PCA embeddings. The centroid distance quantifies how far the flight samples' circadian fingerprints deviate from their ground controls in the reduced space. We then tested whether centroid distance correlated with the absolute phase shift magnitude using Spearman's rank correlation.

#### Gene set enrichment analysis

To connect the observed phase shifts to transcriptional consequences, we performed differential expression and gene set enrichment analysis on the four studies with statistically significant phase shifts (OSD-321, OSD-193, OSD-38, OSD-281). Differential expression was assessed using limma-voom [16], appropriate for the normalized count matrices provided by GeneLab (DESeq2 median-of-ratios). For each study, flight samples were compared to ground controls with the default empirical Bayes moderation. Genes with FDR-adjusted p < 0.05 were considered differentially expressed.

Gene set enrichment analysis was performed using fgsea [17] with the fast preranked GSEA algorithm. Gene sets were constructed from the org.At.tair.db annotation package, combining Gene Ontology Biological Process terms (GO BP) and KEGG pathway annotations. Gene sets containing fewer than 15 or more than 500 genes were filtered out, yielding 816 gene sets. Genes were ranked by the t-statistic from the limma analysis. Pathways with FDR-adjusted p < 0.05 were considered significant. We specifically examined the circadian rhythm pathway (GO:0007623) and regulation of circadian rhythm pathway (GO:0042752) across all four studies. An interactive dashboard enabling exploration of all results is available at the companion website (see Data Availability).

### Fork comparison

Fork A (Barker 2023 accessions, 14 studies) and Fork B (all Arabidopsis spaceflight transcriptomics, 23 studies) were analyzed independently to assess how dataset scope affects conclusions. Fork A represents a curated subset focused on key spaceflight experiments, while Fork B provides a comprehensive analysis of all available data.

## Results

### Dataset overview

We analyzed 23 Arabidopsis thaliana spaceflight transcriptomics datasets from NASA GeneLab with usable processed expression data, comprising 603 samples with ChronoGauge predictions. The datasets spanned 13 ecotypes (dominated by Col-0 with 770 samples), 7 tissue types (root, whole seedling, hypocotyl, cotyledon, shoot, leaf, cell culture), and 12 hardware configurations (including EMCS, Veggie, BRIC-PDFU, ABRS, and petri dishes). The studies included 18 RNA-seq and 9 microarray datasets (some studies had both platforms). Light regimes included continuous light (451 samples), dark conditions (188 samples), and various photoperiods (160 samples), with 370 samples having unspecified light regimes in the ISA-Tab metadata. Of the 23 studies with predictions, 18 had paired flight and ground control samples suitable for phase-shift analysis.

### ChronoGauge validation

On ChronoGauge's held-out RNA-seq test data (Graf et al. dataset, two ecotypes), the ensemble of 100 sub-models achieved a mean absolute error of 44.0 minutes (median error: 21.1 minutes, RMSE: 64.8 minutes). This is consistent with the published ChronoGauge performance and confirms that the model can reliably predict circadian time from single transcriptomic snapshots. The circular variance on validation data averaged 0.18, substantially lower than the GeneLab samples (mean 0.49), reflecting the more controlled conditions of the training data (constant light, synchronized seedlings).

### Per-study phase shifts

Four of 18 studies showed statistically significant phase shifts after Watson-Williams testing (Table 2). The largest shift was observed in OSD-321 (−0.43 h, p = 1.8 × 10⁻⁷, whole seedlings, BRIC-PDFU hardware, dark conditions), representing a 26-minute phase advance in flight samples. OSD-38 showed a −0.52 h shift (p = 0.026, whole seedlings, BRIC-PDFU, dark), OSD-193 showed a −0.24 h shift (p = 0.002, roots, Veggie hardware), and OSD-281 showed a −0.17 h shift (p = 0.038, roots, Veggie hardware). All significant shifts were in the negative direction (phase advance). The remaining 14 studies showed non-significant trends, with 11 of 14 also in the negative direction.

### Meta-analysis

Random-effects meta-analysis across 18 studies revealed a small but statistically significant phase advance under spaceflight (pooled shift = −0.09 h, 95% CI [−0.18, −0.001], p = 0.046, Z = −1.99). Heterogeneity was high (Q = 126.2, p < 10⁻¹⁸, I² = 86.5%), reflecting the diversity of study designs, tissues, hardware, and light regimes. The pooled effect corresponds to approximately a 5.4-minute phase advance, which, while small in absolute terms, is consistent in direction across the majority of studies.

### Stratified analysis

Stratification by tissue type revealed that the phase advance was driven primarily by root studies (pooled shift = −0.17 h, 95% CI [−0.24, −0.09], p < 0.001, I² = 8.9%, k = 4), which showed remarkably low heterogeneity. Whole seedling studies showed a larger but non-significant trend (−0.19 h, p = 0.121, I² = 91.8%), while leaf studies showed a non-significant phase delay (+0.19 h, p = 0.094). Hypocotyl and shoot studies showed no significant effects.

Stratification by light regime revealed that the phase advance was significant only in dark-grown samples (pooled shift = −0.20 h, 95% CI [−0.37, −0.03], p = 0.024, I² = 91.8%, k = 6) and in samples with unknown light regimes (−0.09 h, p = 0.039, I² = 59.6%, k = 7). Continuous light samples showed no significant shift (+0.01 h, p = 0.932), and LD 16:8 samples showed a non-significant trend toward delay (+0.14 h, p = 0.086).

Stratification by hardware was limited by the number of studies per hardware type. The Veggie system (5 studies) showed a non-significant trend (−0.08 h, p = 0.327), while BRIC-PDFU studies showed larger but non-significant shifts due to high heterogeneity.

### Fork A versus Fork B

The forked analysis demonstrated consistent conclusions regardless of dataset scope. Fork A (10 studies in meta-analysis) yielded a pooled phase shift of −0.09 h (95% CI [−0.17, −0.02], p = 0.019, I² = 74.4%). Fork B (18 studies in meta-analysis) yielded a pooled phase shift of −0.09 h (95% CI [−0.18, −0.001], p = 0.047, I² = 86.5%). The point estimates were nearly identical, though Fork A achieved slightly higher statistical significance due to lower heterogeneity from the more curated dataset. This consistency suggests that the circadian phase advance under spaceflight is a robust finding not driven by the inclusion of specific studies.

### Circadian fingerprint

Flight samples showed slightly lower mean circular variance than ground controls (0.488 vs 0.507, difference = −0.019), suggesting marginally more coherent ensemble predictions under spaceflight. This counterintuitive finding may reflect more uniform growth conditions in the controlled spaceflight environment, or it may indicate that the clock output becomes more tightly regulated under stress. However, the difference is small and its biological significance is uncertain.

### Circadian trajectory analysis

Dimensionality reduction on the 100-dimensional ChronoGauge sub-model prediction vectors revealed that the multivariate circadian fingerprint captures phase disruption magnitude. In the t-SNE embedding (Supplementary Figure S7), samples clustered partially by condition and study, with flight and ground centroids separated to varying degrees across studies. The Euclidean distance between flight and ground centroids in t-SNE space correlated strongly with the absolute phase shift magnitude (Spearman ρ = 0.633, p = 0.005, N = 18 studies). A similar correlation was observed in PCA space (ρ = 0.628, p = 0.005). Studies with the largest centroid distances—OSD-321 (12.7), OSD-217 (11.8), and OSD-193 (5.3)—were among those with the largest or most significant phase shifts. This demonstrates that the ensemble's internal prediction structure encodes biologically meaningful information about circadian disruption beyond the aggregate phase estimate.

### Gene set enrichment of phase-shifted studies

To determine whether the observed phase shifts translate into coherent transcriptional changes, we performed limma-voom differential expression and fgsea gene set enrichment on the four studies with significant phase shifts. Differential expression yielded 15,253 DEGs (FDR < 0.05) in OSD-321, 3,481 in OSD-38, 2,603 in OSD-193, and 2,454 in OSD-281, confirming substantial transcriptomic remodeling in these studies.

Gene set enrichment across 816 GO BP and KEGG pathways revealed that OSD-193—the root study with a −0.24 h phase shift (p = 0.002)—showed significant downregulation of the circadian rhythm pathway (GO:0007623, NES = −2.25, FDR = 0.004) (Supplementary Figure S8, Supplementary Table S3). The leading-edge genes driving this enrichment included core clock components LHY (AT1G01060), TOC1 (AT5G61380), PRR5 (AT5G24470), CCA1 (AT2G46830), PRR7 (AT2G46790), and TIC (AT3G46640). This directly connects the ChronoGauge-detected phase advance to reduced expression of the oscillator machinery itself. The other three studies did not show significant circadian pathway enrichment despite having significant phase shifts, suggesting that phase disruption and transcriptional downregulation of clock genes are distinct phenomena that co-occur in some but not all contexts.

## Discussion

This study provides the first systematic quantification of circadian phase disruption in spaceflight-grown plants using deep learning-based circadian time prediction. Our meta-analysis of 18 Arabidopsis spaceflight transcriptomics studies reveals a small but statistically significant circadian phase advance under spaceflight conditions (pooled shift ≈ −0.09 h, approximately 5 minutes). While modest in magnitude, this effect is consistent in direction across the majority of studies and reaches statistical significance in both the curated Fork A and comprehensive Fork B analyses.

The tissue-stratified analysis provides the most informative finding: root studies showed a significant phase advance of −0.17 h (p < 0.001) with remarkably low heterogeneity (I² = 8.9%), suggesting a tissue-specific circadian response to spaceflight. This is consistent with the known sensitivity of root circadian rhythms to mechanical and gravitational signals [14]. Roots are the primary site of gravity sensing in plants, and the removal of gravitational cues in microgravity may directly affect the root circadian oscillator through altered PIN protein cycling and auxin transport rhythms [15].

The light regime stratification revealed that the phase advance was significant only in dark-grown samples (−0.20 h, p = 0.024), while continuous light samples showed no effect. This is biologically plausible: in the absence of light entrainment cues, the circadian clock relies entirely on internal oscillators, which may be more susceptible to disruption by microgravity-induced changes in cellular signaling. Under continuous light, the strong photic entrainment may override any spaceflight-induced phase perturbation. This finding has implications for spaceflight experiment design: studies using dark or dim-light conditions may be more likely to reveal circadian disruption than those using continuous illumination.

The high heterogeneity observed in the overall meta-analysis (I² = 86.5%) reflects the substantial diversity in study designs, including different hardware systems, ecotypes, tissues, light regimes, and growth durations. This heterogeneity is both a limitation and a strength: it demonstrates that the circadian phase advance is not driven by a single experimental condition, but it also means that the pooled effect size should be interpreted as an average across diverse conditions rather than a precise estimate applicable to any specific scenario.

The circadian trajectory analysis provides an independent line of evidence supporting the biological validity of the ChronoGauge predictions. The strong correlation between t-SNE centroid distance and phase shift magnitude (ρ = 0.633, p = 0.005) demonstrates that the ensemble's internal prediction structure—not just the aggregate circular mean—encodes meaningful information about circadian disruption. Studies where flight and ground samples separate widely in the sub-model prediction space are precisely those with the largest phase shifts, suggesting that the 100-dimensional circadian fingerprint captures the multidimensional signature of clock perturbation.

The gene set enrichment analysis further strengthens the connection between predicted phase disruption and transcriptional reality. In OSD-193, the significant downregulation of the circadian rhythm pathway (GO:0007623, NES = −2.25, FDR = 0.004), driven by leading-edge genes including LHY, TOC1, PRR5, CCA1, PRR7, and TIC, directly links the ChronoGauge-detected phase advance to reduced expression of the core oscillator machinery. This is notable because ChronoGauge predicts phase from the full transcriptome, not solely from clock genes; the concordance between predicted phase shift and clock gene pathway downregulation in an independent analysis provides orthogonal validation. However, the absence of circadian pathway enrichment in the other three phase-shifted studies (OSD-321, OSD-38, OSD-281) indicates that phase disruption can occur without coordinated transcriptional suppression of clock genes, possibly reflecting different mechanisms of clock perturbation—for example, altered post-translational regulation or changes in clock output genes rather than the core oscillator.

Several limitations should be noted. First, ChronoGauge was trained on data from plants grown under constant light, while many GeneLab studies used dark or photoperiod conditions. The model's predictions for non-constant-light samples represent extrapolations, and their accuracy may be reduced. The higher circular variance observed in GeneLab samples (mean 0.49) compared to validation data (mean 0.18) likely reflects this domain shift. Second, most GeneLab studies lack circadian time ground truth, preventing direct validation of predictions on spaceflight samples. Third, the cross-sectional nature of the data means we cannot distinguish between true circadian phase shifts and changes in the expression of clock genes that do not reflect actual phase changes. Fourth, the metadata for light regime was incomplete for 370 samples (32%), limiting the power of the light regime stratification. Finally, the small effect sizes (sub-1-hour shifts) are near the limits of ChronoGauge's resolution (MAE = 44 minutes), and individual study results should be interpreted with caution.

Despite these limitations, our approach offers several advantages over traditional circadian analysis methods. By predicting CT from single samples, we can analyze cross-sectional data that would otherwise be uninformative for circadian studies. The ensemble approach provides a built-in uncertainty metric (circular variance) that flags samples where the model is less confident. And the large sample size enabled by analyzing all available GeneLab data provides statistical power that individual studies cannot achieve.

## Conclusions

We applied ChronoGauge, a deep learning circadian time predictor, to 23 Arabidopsis spaceflight transcriptomics datasets from NASA GeneLab. Our random-effects meta-analysis of 18 studies with paired flight and ground controls revealed a small but significant circadian phase advance under spaceflight (−0.09 h, p = 0.046). The effect was strongest in root tissue (−0.17 h, p < 0.001) and in dark-grown samples (−0.20 h, p = 0.024), suggesting that circadian disruption is tissue-specific and modulated by light regime. The forked analysis demonstrated that this finding is robust to dataset scope. Circadian trajectory analysis showed that the multivariate sub-model prediction fingerprint correlates with phase shift magnitude (ρ = 0.633, p = 0.005), and gene set enrichment linked the phase advance in OSD-193 to significant downregulation of the circadian rhythm pathway (GO:0007623, FDR = 0.004) with core clock genes among the leading edge. These results provide the first systematic evidence that spaceflight alters plant circadian clock phase and suggest that future spaceflight experiments should consider light regime and tissue type when assessing circadian disruption. The application of deep learning-based circadian time prediction to retrospective transcriptomics data opens new avenues for extracting circadian information from the wealth of existing cross-sectional datasets in public repositories. An interactive dashboard for exploring all results is available at the companion website.

## References

[1] Covington, M.F., Maloof, J.N., Straume, M., Kay, S.A. & Harmer, S.L. (2008) Genome-wide analysis of the Arabidopsis clock reveals complex regulation. *Genome Biology* 9, R167.

[2] Harmer, S.L., Hogenesch, J.B., Straume, M., Chang, H.S., Han, B., Zhu, T., Wang, X., Kreps, J.A. & Kay, S.A. (2000) Orchestrated transcription of key pathways in Arabidopsis by the circadian clock. *Science* 290, 2110–2113.

[3] Pokhilko, A., Fernández, A.P., Edwards, K.D., Southern, M.M., Halliday, K.J. & Millar, A.J. (2012) The clock gene circuit in Arabidopsis includes a repressilator with additional feedback loops. *Molecular Systems Biology* 8, 574.

[4] Dodd, A.N., Salathia, N., Hall, A., Kévei, E., Tóth, R., Nagy, F., Hibberd, J.M., Millar, A.J. & Webb, A.A.R. (2005) Plant circadian clocks increase photosynthesis, growth, survival, and competitive advantage. *Science* 309, 630–633.

[5] Paul, A.L., Amalfitano, C.E. & Ferl, R.J. (2012) Plant growth strategies are remodeled by spaceflight. *BMC Plant Biology* 12, 232.

[6] Nayak, K.K., et al. (2020) Spaceflight-induced changes in the Arabidopsis circadian clock. [Reference to be verified]

[7] Kiss, J.Z., Kumar, P., Millar, K.D.L., Edelmann, R.E. & Correll, M.J. (2019) Operations of a spaceflight experiment to study plant tropisms. *Gravitational and Space Research* 7, 1–14.

[8] Barker, R., et al. (2023) A meta-analysis of Arabidopsis thaliana spaceflight transcriptome data. *npj Microgravity* 9, [article number].

[9] Hughes, M.E., Hogenesch, J.B. & Kornacker, K. (2010) JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data. *Journal of Biological Rhythms* 25, 372–380.

[10] Thaben, P.F. & Westermark, P.O. (2014) Detecting rhythms in time series with RAIN. *Journal of Biological Rhythms* 29, 391–400.

[11] Reynolds, W., et al. (2025) ChronoGauge: a deep learning approach for circadian time prediction from single transcriptomic samples. *Nature Communications* [details to be verified].

[12] Watson, G.S. & Williams, E.J. (1956) On the construction of significance tests on the circle and the sphere. *Biometrika* 43, 344–352.

[13] Viechtbauer, W. (2010) Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software* 36, 1–48.

[14] Mancuso, S., Barlow, P.W., Volkmann, D. & Baluška, F. (2019) Gravity as a morphogenetic factor in plants. *Journal of Plant Growth Regulation* 38, 1–12.

[15] Rakusová, H., Gallego-Bartolomé, J., Vanstraelen, M., Robert, H.S., Alabadí, D. & Friml, J. (2011) Polarization of PIN3-dependent auxin transport for hypocotyl gravitropic response in Arabidopsis thaliana. *The Plant Journal* 67, 817–826.

[16] Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W. & Smyth, G.K. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research* 43, e47.

[17] Korotkevich, G., Sukhov, V., Budin, N., Shpak, B., Artyomov, M.N. & Sergushichev, A. (2021) Fast gene set enrichment analysis. *bioRxiv* 2021.06.15.448424.

## Data availability

All data are publicly available from NASA GeneLab (https://genelab.nasa.gov). ChronoGauge pre-trained models are available from HuggingFace. Analysis code is available at [GitHub URL] and archived at [Zenodo DOI]. An interactive dashboard for exploring all results—including per-sample predictions, forest plots, circadian fingerprints, clock gene heatmaps, and t-SNE trajectories—is deployed at [GitHub Pages URL].

## Code availability

All analysis code is available at [GitHub URL] under the MIT license. A Docker container for reproducible analysis is available at [DockerHub URL].

## Acknowledgements

This research used data from NASA's GeneLab database. We thank the GeneLab team and all data contributors. We thank the ChronoGauge developers for making pre-trained models publicly available.

## Author contributions

R.B. conceived the study, designed the analysis, performed all computational work, and wrote the manuscript.

## Competing interests

The authors declare no competing interests.
