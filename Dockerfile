FROM eddelbuettel/r2u:24.04
LABEL org.name="Alamar Biosciences, USA"

ENV HOME="/workingDir"
ARG GH_PAT
ARG SOURCE_BRANCH
ARG NULISASEQAQ_BRANCH

WORKDIR ${HOME}

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends pandoc chromium-browser libssl-dev pngquant && \
    rm -rf /var/lib/apt/lists/*

RUN install.r future pagedown uuid XML xml2 ggrepel emmeans fields optparse pheatmap rcompanion \
        corrplot pander gdata dplyr testthat qs lmerTest logger devtools remotes knitr tidyverse \
        ggplot2 htmltools DT plotly lubridate scales httr stringr ggpubr openxlsx renv coin \
        fs reactable janitor future.apply future.callr showtext svglite patchwork

# Install ggalt from RSPM snapshot (removed from CRAN on 2025-08-02)
RUN Rscript -e "install.packages('ggalt', repos = 'https://packagemanager.posit.co/cran/2025-08-01')"

RUN Rscript -e "remotes::install_version('kableExtra', version = '1.3.4', repos = 'https://cran.rstudio.com/')" && \
    Rscript -e "devtools::install_github('Alamar-Biosciences/PCAtools', auth_token = '${GH_PAT}')" && \
    Rscript -e "library(methods); if (!requireNamespace('PCAtools', quietly = TRUE)) stop('Package PCAtools not installed')" && \
    Rscript -e "library(methods); if (!requireNamespace('kableExtra', quietly = TRUE)) stop('Package kableExtra not installed')" && \
    Rscript -e "remotes::install_github('Alamar-Biosciences/ComplexHeatmap', ref='master', auth_token = '${GH_PAT}')" && \
    Rscript -e "library(methods); if (!requireNamespace('ComplexHeatmap', quietly = TRUE)) stop('Package ComplexHeatmap not installed')" && \
    Rscript -e "devtools::install_github('Alamar-Biosciences/NULISAseqAQ', ref='${NULISASEQAQ_BRANCH}', auth_token = '${GH_PAT}')" && \
    Rscript -e "library(methods); if (!requireNamespace('NULISAseqAQ', quietly = TRUE)) stop('Package NULISAseqAQ not installed')"

RUN mkdir -p ${HOME}/NULISAseqR
COPY . ${HOME}/NULISAseqR/.

RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('/workingDir/NULISAseqR', repos=NULL, type='source')" && \
    Rscript -e "library(methods); if (!requireNamespace('NULISAseqR', quietly = TRUE)) stop('Package NULISAseqR not installed')"

EXPOSE 8000
ENV PATH "${HOME}:$PATH"

CMD ["bash"]
