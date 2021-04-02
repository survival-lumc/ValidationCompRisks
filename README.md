# External validation of the performance of competing risks prediction models: a guide through modern methods

R Code repository for the manuscript 'External validation of the performance of competing risks prediction models: a guide through modern methods' (in preparation).

The repository contains the following code:

+ [Prediction_FG](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_FG.md) illustrates how to develop and validate a competing risks prediction model using the Fine and Gray subdistribution hazard regression model. The RMarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_FG.Rmd).  

+ [Prediction_CSC](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC.md) illustrates how to develop and validate a competing risks prediction model using cause-specific hazards regression models.he RMarkdown source code (.Rmd) is [here](https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC.Rmd)

Data are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/Data).  
Additional functions useful to develop and validate competing risks predition models are available [here](https://github.com/survival-lumc/ValidationCompRisks/tree/main/Functions).

## Usage

You can either download a zip file containing the directory, or you can clone it by using

```bash
git clone https://github.com/survival-lumc/ValidationCompRisks.git
```

In either case, you can then use the `ValidationCompRisks.Rproj` file to open
and Rstudio session in the directory you have just downloaded. You may then knit
both rmarkdown files, or run them line-by-line.

## Contributions

| Name                                                         | Affiliation                           | Role              |
| ------------------------------------------------------------ | ------------------------------------- | ----------------- |
| [Daniele Giardiello](https://github.com/danielegiardiello/)  | The Netherlands Cancer Institute (NL) | Author/maintainer |
| [Edouard Bonneville](https://www.lumc.nl/org/bds/medewerkers/1968807) | Leiden University Medical Center (NL) | Code reviewer     |