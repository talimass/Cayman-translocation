## Zooxantellae identification

### Zooxantellae proteomes
For the analysis, 12 proteomes of zooxantellae species were downloaded from open databases: 

|   prefix and link     |   symbiont species                   | former clade according to [LaJeunesse et al., 2018](https://www.cell.com/current-biology/fulltext/S0960-9822(18)30907-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982218309072%3Fshowall%3Dtrue)|
|---------------|--------------------------------------|---|
| [Dtrenchii](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_963970005.1/)     |  Durusdinium trenchii                       | former clade D |
| [symb](https://marinegenomics.oist.jp/symb/viewer/info?project_id=21)          | Breviolum minutum                    | former clade B |
| [syma](https://marinegenomics.oist.jp/symb/viewer/download?project_id=37)          | Symbiodinium clade A3                |   |
| [Stri](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium tridacnidorum           | former clade A  |
| [SPilosum](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905231905.1/)      | Symbiodinium pilosum                 | former clade A  |
| [Snec](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium necroappetens           | former clade A  |
| [Snat](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium natans                  | former clade A  |
| [SMicReefgen](http://smic.reefgenomics.org/)   | Symbiodinium microadriaticum     | former clade A  |
| [Slin](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium linuacheae              | former clade A  |
| [SKawaF](http://symbs.reefgenomics.org/download/)        | Fugacium kawagutii                   | former clade F  |
| [SGoreauC](http://symbs.reefgenomics.org/download/)      | Cladocopium goreau                   | former clade C1  |
| [CladocPluteaC](https://marinegenomics.oist.jp/symb/viewer/download?project_id=40) | Cladcopium C15                       |   |

### DIAMOND run

The prefixes were added to the headers of each fasta file with `sed -i "s/^>/>prefix_/" file.fa`. Fasta files were concatenated with `cat`, and this file was used for a DIAMOND database construction.

Since DIAMOND can't handle paired reads, unmapped R1 and R2 holobiont reads were concatenated into one fastq file per sample with `cat` command before the analysis. Filenames were shortened with `for f in *_combined.fastq*; do mv -- "$f" "${f/_combined.fastq/}"; done` command. 

### Further analysis

After the mapping with DIAMOND, the number of reads aligned to each symbiont proteome in the resulting *.m8 files was counted and visualized with R scripts. R scripts were based on [the script](https://github.com/fscucchia/Pastreoides_development_depth/blob/main/Species_Identification/diamond_heatmap.R) of [F. Scucchia](https://github.com/fscucchia).



