FROM 033947443408.dkr.ecr.us-west-1.amazonaws.com/nulisa-base:v1.0

EXPOSE 8000

WORKDIR /workingDir
RUN mkdir -p /workingDir/NULISAseqR
COPY . /workingDir/NULISAseqR/.

RUN chmod 755 /workingDir/NULISAseqR/R/*.R
RUN apk update && \
    apk --no-cache --update-cache add openssl-dev make libc-dev libxml2 libxml2-dev libxext libxt g++ && \
    rm -rf /var/cache/apk/*
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('future')";
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('pagedown')";
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('uuid')";
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('XML')";
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('xml2')";
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('fields')"; 
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages('/workingDir/NULISAseqR', repos=NULL, type='source')";
ENV PATH "/workingDir:$PATH"
COPY main.R /workingDir/.

CMD ["bash"]
