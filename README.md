# Aritcleone: Eicosanoids as Mediators of Hypertension and Hypertension-Related Cardiovascular Disease

This repository contains computational notebook for the first scientific article of my thesis. 

File                        | Purpose
--------------------------- | -----------------------------------
articleone.Rmd              | Generates docx of the Results
articleone-supplement.Rmd   | Generates docx for the Supplement
articleone-presentation.Rmd | Generates pptx for the presentation
plotmziddistributions.R     | Draws some distribution plots


## Dependencies

R Markdown files benefit from some general purpose functions of *TurkuMetabolite*.

TurkuMetabolite can be installed in Rstudio by SSH

```
cred <- git2r::cred_ssh_key("~/.ssh/id_rsa.pub", "~/.ssh/id_rsa")                           
devtools::install_git(repo, credentials = cred) 
```

Mostly for dataimport some functions of *biodataCore* are also used. Namely *TurkuMetabolite* has dependencies for *biodatacore/biodatacoreImport*, *biodatacore/biodatacoreUtils*, and *biodatacore/biodatacoreSmallMolecules*.



