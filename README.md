# Validation of prediction models in presence of competing risks: a guide through modern methods

R Code repository for the manuscript ['Validation of prediction models in presence of competing risks: a guide through modern methods'](https://www.bmj.com/content/377/bmj-2021-069249) published in BMJ.

The repository contains the following code:


+ **[Prediction_CSC_minimal.R](Prediction_CSC_minimal.R) : the companion R script for the manuscript. This script reproduces all main tables and figures of the manuscript. The file evaluates the performance of a cause specific hazards prediction model. Data are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/Data).**

+ [Prediction_CSC.md](Prediction_CSC.md) : a markdown document containing a more in-depth version of the script, with details on model development, descriptive tables and supplementary plots. This document requires installing several add-on packages that are not needed to run the minimal script. The RMarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC.Rmd). Additional functions useful for this script are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/R). 


+	Additional code to alternatively develop a competing risk prediction model using the subdistribution hazard approach (Fine & Gray) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH.md). The Rmarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH.Rmd). A more concise R source code (.R) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Development_SDH_minimal.R).

+ [sharing_CSC_model.R](sharing_CSC_model.R) : example/template of how to share a cause-specific hazards prediction model for external validation, without having to share the original development data.


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
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | The Netherlands Cancer Institute (NL) <br /> Leiden University Medical Center (NL) <br /> EURAC research (IT) | Author .Rmd files/maintainer                    |
| [Edouard Bonneville](https://www.linkedin.com/in/edbonneville/?originalSubdomain=nl) | Leiden University Medical Center (NL) | Author minimal .R and CSC sharing scripts/review of .Rmd scripts |
| [Nan van Geloven](https://www.universiteitleiden.nl/medewerkers/nan-van-geloven#tab-1) | Leiden University Medical Center (NL) | Review of both .R and .Rmd scripts              |
| [Maarten van Smeden](https://www.umcutrecht.nl/en/research/researchers/van-smeden-maarten-m) | University Medical Centre Utrecht (NL) |Review of .Rmd script     |

