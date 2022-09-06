# NULISAseqR
NULISAseq R package

## How to install

1. Create a personal access token at 
[https://github.com/settings/tokens](https://github.com/settings/tokens)

2. Run the following R code, where `tokenstring` is the personal access token:
```
    install.packages('devtools')
    devtools::install_github('Alamar-Biosciences/NULISAseqR',
                              ref = 'main',
                              auth_token = 'tokenstring'
                              )
```

3. Load the package in R with `library(NULISAseqR)`.
