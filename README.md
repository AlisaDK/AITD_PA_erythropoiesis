# AITD_PA_erythropoiesis: Thyroid function, pernicious anemia, and erythropoiesis: a two-sample Mendelian randomization study

This repository contains data and code used for "Thyroid function, pernicious anemia, and erythropoiesis: a two-sample Mendelian randomization study" (doi: 10.1093/hmg/ddac052)

# ABSTRACT
Autoimmune thyroid disease (AITD) and pernicious anemia (PA) often coexist, but the directionality is unknown. In a two-sample Mendelian randomization (MR) analysis, using summary statistics from large genome-wide association studies in Europeans (N = 49 269-755 406), we examined the genetic associations between thyroid function, pernicious anemia and markers of erythropoiesis. We performed inverse variance weighted random-effects MR, several sensitivity MR analyses, and bidirectional MR and MR Steiger for directionality. AITD and PA were associated bidirectionally (P ≤ 8E-6). Neither euthyroid thyroid stimulating hormone (TSH) nor free thyroxine (FT4) were causally associated with PA. One standard deviation increase in euthyroid FT4 regulated by genetic variants in deiodinases 1 and 2 genes (DIO1/DIO2), corresponding to low-normal free triiodothyronine (FT3) levels, was causally associated with a pernicious/macrocytic anemia pattern, i.e. decreased erythrocyte counts (rank-based inverse normal transformed β = -0,064 [95% confidence interval: -0,085,-0,044], P = 8E-10) and hemoglobin (-0.028 [-0.051,-0.005], P = 0.02) and increased mean corpuscular hemoglobin (0.058 [0.025,0.091], P = 5E-4) and mean corpuscular volume levels (0.075 [0.052,0.098], P = 1E-8). Meanwhile, subclinical hyperthyroidism mirrored that pattern. AITD was causally associated with increased erythrocyte distribution width (P = 0.007) and decreased reticulocyte counts (P ≤ 0.02), whereas high-normal FT4 regulated by DIO1/DIO2 variants was causally associated with decreased bilirubin (-0.039 (-0.064,-0.013), P = 0.003). In conclusion, the bidirectional association between AITD and PA suggests a shared heritability for these two autoimmune diseases. AITD was causally associated with impaired erythropoiesis and not autoimmune hemolysis. Additionally, in euthyroid individuals, local regulation of thyroid hormones by deiodinases likely plays a role in erythropoiesis.


# Brief description of folder and file contents

The following folders contain:

- `data-raw/`: summary statistics data input
- `data-raw/b12_meta_ukbb_estbb_finngen_061020.out`: please download from from http://www.geenivaramu.ee/tools/pernicious_anemia_Laisketal2021_sumstats.gz (GWAS of 2,166 pernicious anemia cases and 659,516 controls from population-based biobanks: Estonian Biobank, UK Biobank and FinnGen study)
- `data-raw/Autoimmune_thyroid_disease_Meta_results_200313.txt`: please download from https://www.nature.com/articles/s41586-020-2436-0, https://www.decode.com/summarydata/ (GWAS of 30,234 AITD cases and 755,172 controls from Iceland and UK Biobank)
- `data/`: generated data
- `R/`: R code for analyses