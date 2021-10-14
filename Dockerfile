FROM quay.io/aptible/ubuntu:18.04
RUN apt-get update && \
    apt-get install -y python3.8 python3-pip python3.8-dev wget curl
RUN /usr/bin/env python3 -m pip install --upgrade pip && \
    /usr/bin/env python3 -m pip install \
    dash-bootstrap-components==0.9.2 \
    dash-core-components==1.9.1 \
    dash-extensions \
    dash-html-components==1.0.3 \
    dash-renderer==1.4.0 \
    dash-table==4.6.2 \
    dash==1.11.0 \
    direct-redis \
    matplotlib \
    numpy==1.19.5 \
    pandas \
    plotly \
    plotnine>=0.8.0 \
    psutil \
    seaborn==0.11.1 && \
    /usr/bin/env python3 -m pip install \
    scikit-bio==0.5.6 \
    scikit-learn==0.24.1 \
    scikit-image \
    biopython \
    watchdog \
    xarray==0.16.2
RUN mkdir /usr/local/redis && \
    cd /usr/local/redis && \
    wget http://download.redis.io/redis-stable.tar.gz && \
    tar xvzf redis-stable.tar.gz && \
    cd redis-stable && \
    make && \
    make install
RUN /usr/bin/env python3 -m pip install \
    git+https://github.com/FredHutch/menu-driven-figure.git@e29bef4
ADD app/ /usr/bin/local/
RUN ln -s /usr/bin/local/gig-map /usr/local/bin/ && \
    gig-map --help \
    gig-map-cli --help
RUN mkdir /work
WORKDIR /work