FROM 033947443408.dkr.ecr.us-west-1.amazonaws.com/nulisa-base:v1.0

EXPOSE 8000

WORKDIR /workingDir
RUN mkdir -p /workingDir/NULISAseqR
COPY . /workingDir/NULISAseqR/.

RUN chmod 755 /workingDir/NULISAseqR/R/*.R
RUN apk add --no-cache make libc-dev libxml2-dev libxext libxt g++
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages(c('future','uuid','XML','fields', 'xml2')); install.packages('/workingDir/NULISAseqR', repos=NULL, type='source')"
ENV PATH "/workingDir:$PATH"
COPY main.R /workingDir/.

CMD ["bash"]
