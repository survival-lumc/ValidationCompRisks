# External validation of the performance of competing risks prediction models: a guide through modern methods

R Code repository for the manuscript 'External validation of the performance of competing risks prediction models: a guide through modern methods' (in preparation).

The repository contains the following code:

+ [Prediction_CSC_minimal.R](Prediction_CSC_minimal.R) : the companion (minimal) script for the manuscript, illustrating external validation of a prediction model. The file uses a cause specific hazards prediction model. To reproduce all mean tables and figures of the manuscript, this script is sufficient. 

+ [Prediction_CSC.md](Prediction_CSC.md) : a markdown document containing a more in-depth version script, with details on model development, descriptive tables and plots. The RMarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC.Rmd).

+ Additional code to alternatively develop a competing risk prediction model using the subdistribution hazard approach (Fine & Gray) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH.md). The Rmarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH.Rmd). A more concise R source code (.R) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH_minimal.R).

+ [sharing_CSC_model.R](sharing_CSC_model.R) : example/template of how to share a cause-specific hazards prediction model for external validation, without having to share the original development data.


Data are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/Data).  Additional functions useful to develop and validate competing risks prediction models are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/R).

## Usage

If you are git user, you can clone the directory by using

```bash
git clone https://github.com/survival-lumc/ValidationCompRisks.git
```

Otherwise, you can simply download a zip file containing the directory by clicking Code -> Download ZIP at the top-right of this Github page. Extract the zipped files to a directory of your choice.

Afterwards, you can double-click the `ValidationCompRisks.Rproj` file to open an Rstudio session in the directory you have just downloaded. This will ensure all file-paths called in the files are maintained. The minimal script and the .Rmd files can now be executed.

## Contributions

| Name                                                         | Affiliation                           | Role                                            |
| ------------------------------------------------------------ | ------------------------------------- | ----------------------------------------------- |
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | The Netherlands Cancer Institute (NL) | Author .Rmd files/maintainer                    |
| [Edouard Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807) | Leiden University Medical Center (NL) | Author minimal .R and CSC sharing scripts/review of .Rmd scripts |
| [Nan van Geloven](https://www.lumc.nl/org/bds/medewerkers/1216536?setlanguage=English&setcountry=en) | Leiden University Medical Center (NL) | Review of both .R and .Rmd scripts              |

