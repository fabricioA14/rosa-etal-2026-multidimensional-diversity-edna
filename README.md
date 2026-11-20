<div align="center">

# Rethinking stream fish assessments: environmental DNA broadens the multidimensional diversity of assemblages in Eastern Amazonia

**Fabricio dos Anjos Santa Rosa**<sup>1,2&#42;</sup>, **Derlan José Ferreira Silva**<sup>1&#42;</sup>, Cintia Oliveira Carvalho<sup>1,3</sup>, Silvia Britto Barreto<sup>2</sup>, Silvana Melo Sviggum<sup>4</sup>, Anna Karolina Oliveira de Queiroz<sup>1,2</sup>, Enki Gautier<sup>1</sup>, Marcia da Silva Anjos<sup>1</sup>, Jarl Andreas Anmarkrud<sup>5</sup>, João Bráullio de Luna Sales<sup>1</sup>, Quentin Mauvisseau<sup>5</sup>, Hugo de Boer<sup>5</sup>, Jonathan Stuart Ready<sup>1</sup>

</div>

<sup>1</sup>Group for Integrated Biological Investigation, Center for Advanced Studies of Biodiversity, Federal University of Pará. Av. Perimetral da Ciência, km 01, PCT‐Guamá, Lot 11, Belém-PA. ZIP 66075‐750, Brazil<br>
<sup>2</sup>Vale Institute of Technology Sustainable Development (ITV DS), Belém-PA, Brazil<br>
<sup>3</sup>The School of Business, Innovation and Sustainability (FIH), Halmstad University, Halmstad, Sweden<br>
<sup>4</sup>Oslo University Hospital, Oslo, Norway<br>
<sup>5</sup>Natural History Museum of Oslo, University of Oslo, Norway

<sup>&#42;</sup>Joint first authors

Corresponding author: Fabricio dos Anjos Santa Rosa (fabriciodosanjos14@gmail.com)

## Study

Fish assemblages of 12 first- to third-order streams of the Potiritá River sub-basin (Capim River basin, eastern
Amazonia) were sampled with seine nets and hand sieves and with water eDNA (COI metabarcoding; nine 400 mL water
samples per stream) in four campaigns: August 2022 and August 2023 (dry season), February 2023 and February 2024
(rainy season). For each stream–campaign unit, the taxonomic (richness), phylogenetic (Faith's PD) and functional
(FRic) diversity recovered by each method, alone and combined, was compared within and between seasons.

## Contents

| File | Description |
|---|---|
| `01_phylogeny_pool_tree.R` | Step 1. Pruned time-calibrated tree of the regional pool (119 taxa) and taxon–tip correspondence (Table S1). |
| `02_phylogenetic_diversity_PD.R` | Step 2. Faith's PD per stream–campaign unit and sampling category. |
| `03_functional_richness_FRic.R` | Step 3. Trait transformation and correlation filter (Table S2), functional space (PCoA) and FRic. |
| `04_indices_per_unit.R` | Step 4. Richness, PD and FRic standardized by the regional pool. |
| `05_complementarity_FigS1_S2_S3.R` | Step 5. Exclusive contribution of each method; Figures S1–S3. |
| `06_season_analyses_Fig2_Fig3_FigS4.R` | Step 6. Analyses within and between seasons; Figures 2, 3 and S4; Table S3. |

The eDNA occurrences (82 taxa) come from the read processing and curation pipeline in
[rosa-etal-2026-edna-vs-traditional-amazon-streams](https://github.com/fabricioA14/rosa-etal-2026-edna-vs-traditional-amazon-streams).

## How to run

This repository holds only the scripts. The data they use will be made available upon acceptance of the associated
article. To run the scripts, place these files in this folder:

- `SuppInfo_S1.xlsx` (Appendix S1 of the article; occurrences in sheet `3_analysis_input`);
- `pool_taxa_traits.xlsx` (taxon list, tree tips and the 33 morphological traits);
- `South-American-freshwater-fish-Fishes.zip` (time-calibrated phylogeny of Cassemiro et al., 2023,
  https://github.com/fcassemiro/South-American-freshwater-fish);
- `silhouettes/` (fish outlines drawn in Figures S2 and S3).

Then run each script from this folder, in order:

```
Rscript 01_phylogeny_pool_tree.R
Rscript 02_phylogenetic_diversity_PD.R
Rscript 03_functional_richness_FRic.R raw
Rscript 04_indices_per_unit.R
Rscript 05_complementarity_FigS1_S2_S3.R
Rscript 06_season_analyses_Fig2_Fig3_FigS4.R
```

Outputs are written to `phylogeny/`, `functional/raw/` and `integrated/`. Section (G) of Table S3 (sensitivity analysis with
untransformed traits) is obtained by running steps 3–6 in a copy of this folder with `NO_LOG=1` set for step 3
(`NO_LOG=1 Rscript 03_functional_richness_FRic.R raw`).

**Software:** R 4.6.1; ape 5.8.1, picante 1.8.4, FD 1.0-12.5, geometry 0.5.2, ggplot2, patchwork, dplyr, tidyr,
openxlsx and readxl. Null models use a fixed seed, so the statistics reported in the manuscript are reproduced exactly.

## Data availability

The dataset is not included in this repository and will be made available upon acceptance of the associated article:
raw sequencing data on Zenodo (https://doi.org/10.5281/zenodo.20977161; under embargo until acceptance) and processed
data (occurrences, read counts, curation decisions and Tables S1–S3) as Appendix S1 of the article (`SuppInfo_S1.xlsx`).

## License

The scripts are released under the MIT License (also in `LICENSE`). The license covers only the scripts; the data are
not part of this repository (see Data availability).

```
MIT License

Copyright (c) 2026 Rosa et al.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

"Rosa et al." refers to all authors listed above.

## Citation

Please cite the article (reference to be added upon publication) and the Zenodo dataset.
