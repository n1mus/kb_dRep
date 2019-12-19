FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.



WORKDIR /kb/module


RUN apt-get update
RUN pip install --upgrade pip==19.3.1


RUN pip install drep==2.3.2


RUN curl --location https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar > mash.tar && \
    tar xf mash.tar && \ 
    mv mash-Linux64-v2.2/mash /usr/local/bin/ && \
    rm -r mash*


RUN apt-get install --yes gcc=4:6.3.0-4 && \
    apt-get install --yes --reinstall zlibc=0.9k-4.3 zlib1g=1:1.2.8.dfsg-5 zlib1g-dev=1:1.2.8.dfsg-5

# TODO version
# dRep says v2.6.3 confirmed works, on its wiki
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


RUN apt-get install --yes hmmer=3.1b2+dfsg-5


RUN apt-get install --yes libbz2-dev=1.0.6-8.1 liblzma-dev=5.2.2-1.2+b1

RUN pip install checkm-genome==1.1.1

# just to get dataRoot pointing to /data/CHECKM_DATA without a 'Permission denied' error
# does this work when re-caching Dockerfile without running refdata scripts?
RUN checkm data setRoot /data/CHECKM_DATA && \
    echo "ls -a /data/CHECKM_DATA" && ls -a /data/CHECKM_DATA && \
    echo "cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG" && cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG


RUN apt-get install g++=4:6.3.0-4 csh=20110502-2.2+b1 --yes

RUN curl --location https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download > mummer3.23.tar.gz && \
tar -vxzf mummer3.23.tar.gz  && \
mv MUMmer3.23 /usr/local/bin && \
rm mummer3.23.tar.gz && \
cd /usr/local/bin/MUMmer3.23/ && \
make check && \
make install

ENV PATH="${PATH}:/usr/local/bin/MUMmer3.23"


RUN curl --location https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz > ANIcalculator_v1.tgz && \
tar vxzf ANIcalculator_v1.tgz && \
rm ANIcalculator_v1.tgz && \
mv ANIcalculator_v1 /usr/local/bin

ENV PATH="${PATH}:/usr/local/bin/ANIcalculator_v1"


# Utilities for manual inspection of Docker container
RUN apt-get install --yes vim tree

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
