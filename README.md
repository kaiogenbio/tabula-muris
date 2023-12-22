# Tabula Muris

## About

The _Tabula Muris_ comprises single-cell transcriptomic data from **100,605 cells** isolated from 20 organs from three female and four male, **C57BL/6JN**, three-month-old mice (10–15 weeks), analogous to 20-year-old humans.

Twenty organs:
1. Aorta
2. Bladder
3. Bone marrow
4. `Brain_Myeloid` (Cerebellum, Cortex, Hippocampus, Striatum)
5. `Brain_Non-Myeloid` (Cerebellum, Cortex, Hippocampus, Striatum)
6. Diaphragm
7. Fat: Brown, Gonadal, Mesenteric, Subcutaneous
8. Heart
9. Kidney
10. Large intestine
11. Limb muscle
12. Liver
13. Lung
14. Mammary gland
15. Pancreas
16. Skin
17. Spleen
18. Thymus
19. Tongue
20. Trachea

## General properties of the data
- As expected, cells from different organs often mixed, with 25 of 54 t-SNE clusters containing at least five cells from distinct organs.
- 1,016 transcription factors expressed in our dataset.

### Method details
- Sequences from the NovaSeq were de-multiplexed using **`bcl2fastq` version 2.19.0.316**.
- Reads were aligned using to the mm10plus genome using **STAR version 2.5.2b** with parameters TK.
- Gene counts were produced using **HTSEQ version 0.6.1p1** with default parameters, except `stranded` was set to `false`, and `mode` was set to `intersection-nonempty`.
- Sequences from the microfluidic droplet platform were de-multiplexed and aligned using **CellRanger version 2.0.1**, available from 10x Genomics with default parameters.

## `raw/`

- `tm_10x_v2.zip` (Tabula Muris sequenced using the 10X Genomics Platform; version 2)
  - source: [figshare:5968960](https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960)
- `tm_ss2.zip` (Tabula Muris sequenced using Smartseq2; version 2)
  - source: [figshare:5829687](https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_Smart-seq2_sequencing_of_FACS_sorted_cells_v2_/5829687)

## `build/`

### `preprocess/ss2`
- Smartseq2 data needs to be merged with annotations file.

### `preprocess/10x`

## Tests

