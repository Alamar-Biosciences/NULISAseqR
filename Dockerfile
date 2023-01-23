FROM 033947443408.dkr.ecr.us-west-1.amazonaws.com/nulisa-base:v1.0

EXPOSE 8000

WORKDIR /workingDir
RUN mkdir -p /workingDir/NULISAseqR
COPY . /workingDir/NULISAseqR/.

RUN chmod 755 /workingDir/NULISAseqR/R/*.R
RUN apk add --no-cache make libc-dev libxml2-dev libxext libxt
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages(c('fields', 'xml2', 'zeallot'))"
RUN R --no-save -e "install.packages('/workingDir/', repos=NULL, type='source')"
ENV PATH "/workingDir:$PATH"

CMD ["bash"]
