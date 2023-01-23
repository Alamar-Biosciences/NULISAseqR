FROM 033947443408.dkr.ecr.us-west-1.amazonaws.com/nulisa-base:v1.0

EXPOSE 8000

WORKDIR /workingDir

COPY ./.git/refs/heads/main ./
COPY ./R/*.R /workingDir/.
COPY main.R ./

RUN chmod 755 ./*.R
RUN R --no-save -e "options(repos=structure(c(CRAN='https://cran.wustl.edu/'))); install.packages(c('fields', 'xml2', 'zeallot'))"
ENV PATH "/workingDir:$PATH"

CMD ["bash"]
