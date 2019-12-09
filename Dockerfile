FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.



WORKDIR /kb/module


RUN apt-get update
RUN pip install --upgrade pip


RUN pip install drep


RUN curl --location https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar > mash.tar && \
    tar xf mash.tar && \ 
    mv mash-Linux64-v2.2/mash /usr/local/bin/ && \
    rm -r mash*


RUN apt-get install --yes gcc && \
    apt-get install --yes --reinstall zlibc zlib1g zlib1g-dev

RUN git clone https://github.com/hyattpd/Prodigal && \
    cd Prodigal/ && \
    make install && \
    cd .. && \
    rm -rf Prodigal


RUN curl --location https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip > pplacer.zip && \
    unzip pplacer.zip && \
    cd pplacer-Linux-v1.1.alpha17/ && \
    mv * /usr/local/bin && \
    cd .. && \
    rm -r pplacer*


RUN apt-get install --yes hmmer


RUN apt-get install --yes libbz2-dev liblzma-dev

RUN pip install checkm-genome && \
    checkm data setRoot /data/CHECKM_DATA


RUN apt-get install --yes vim
RUN apt-get install --yes tree

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
