# NULISAseqR

![AWS CodeBuild](https://codebuild.us-west-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoib045RnFOTFB4Wmo5OHBDTGJySnNJK3dtN2I3a0MwQm96UVZyMnp1anl3cGZtMWs5dVowMVl5TVlLUEw4RnNiZWlscnNTdE5KV2xQSlVyN3YrZUVvYTZRPSIsIml2UGFyYW1ldGVyU3BlYyI6InNtclNBUGloQjJEdytnMUQiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=main)

NULISAseq R package

## How to install

1. Create a personal access token at 
[https://github.com/settings/tokens](https://github.com/settings/tokens). Check the "repo" box.

2. Run the following R code, where `tokenstring` is the personal access token:
```
    install.packages('devtools')
    devtools::install_github('Alamar-Biosciences/NULISAseqR',
                              ref = 'main',
                              auth_token = 'tokenstring'
                              )
```

3. Load the package in R with `library(NULISAseqR)`.
