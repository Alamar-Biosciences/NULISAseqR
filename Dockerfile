FROM eddelbuettel/r2u:22.04
LABEL org.name="Alamar Biosciences, USA"

ENV HOME="/workingDir"

WORKDIR ${HOME}
RUN mkdir -p ${HOME}/NULISAseqR
COPY . ${HOME}/NULISAseqR/.
COPY main.R ${HOME}/.

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends pandoc chromium-browser && \
    rm -rf /var/lib/apt/lists/*

RUN set -eu; \
    install.r future pagedown uuid XML xml2 ggrepel fields plumber optparse kableExtra pheatmap corrplot pander gdata dplyr testthat \
    > /tmp/install.log 2>&1; \
    cat /tmp/install.log; \
    if grep -qi "not available for this version of R" /tmp/install.log; then \
        package=$(grep -Eo '‘[[:alnum:]]+’' /tmp/install.log | tr -d '‘’'); \
        echo "Problem with installing R packages: ${package}"; \
        exit 1; \
    fi; \
    rm -rf /tmp/install.log;

RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('/workingDir/NULISAseqR', repos=NULL, type='source')" && \
    Rscript -e "library(methods); if (!requireNamespace('NULISAseqR', quietly = TRUE)) stop('Package NULISAseqR not installed')"

EXPOSE 8000
ENV PATH "${HOME}:$PATH"

CMD ["bash"]
