## Zooxantellae identification

### Zooxantellae proteomes
For the analysis, 15 proteomes of zooxantellae species were downloaded from open databases: 

|   prefix and link     |   symbiont species                   | comment |
|---------------|--------------------------------------|---|
| [symbd](https://marinegenomics.oist.jp/symbd/viewer/download?project_id=102)         | Durusdinium sp.                      | former clade D |
| [symb](https://marinegenomics.oist.jp/symb/viewer/info?project_id=21)          | Breviolum minutum                    | former clade B |
| [syma](https://marinegenomics.oist.jp/symb/viewer/download?project_id=37)          | Symbiodinium clade A3                |   |
| [Stri](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium tridacnidorum           | former clade A  |
| [SPilosum](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905231905.1/)      | Symbiodinium pilosum                 | former clade A  |
| [Snec](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium necroappetens           | former clade A  |
| [Snat](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium natans                  | former clade A  |
| [SMicUQSCI](https://espace.library.uq.edu.au/view/UQ:f1b3a11)     | Symbiodinium microadriaticum UQ      | former clade A  |
| [SMicUQCass](https://espace.library.uq.edu.au/view/UQ:f1b3a11)    | Symbiodinium microadriaticum UQ Cass | former clade A  |
| [SMicReefgen](http://smic.reefgenomics.org/)   | Symbiodinium microadriaticum RG      | former clade A  |
| [Slin](https://espace.library.uq.edu.au/view/UQ:f1b3a11)          | Symbiodinium linuacheae              | former clade A  |
| [SKawaF](http://symbs.reefgenomics.org/download/)        | Fugacium kawagutii                   | former clade F  |
| [SGoreauC](http://symbs.reefgenomics.org/download/)      | Cladocopium goreau                   | former clade C1  |
| [CladocPluteaC](https://marinegenomics.oist.jp/symb/viewer/download?project_id=40) | Cladcopium C15                       |   |

The prefixes were added to the headers of each fasta file with `sed -i "s/^>/>prefix_/" file.fa`. Fasta files were concatenated with `cat`, and this file was used for a DIAMOND database construction.  

