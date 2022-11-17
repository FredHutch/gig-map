FROM quay.io/aptible/ubuntu:18.04
RUN apt-get update && \
    apt-get install -y python3.8 python3-pip python3.8-dev wget curl jq
RUN python3.8 -m pip install --upgrade pip && \
    python3.8 -m pip install \
    click \
    dash-bootstrap-components==0.9.2 \
    dash-core-components==1.9.1 \
    dash-extensions \
    dash-html-components==1.0.3 \
    dash-renderer==1.4.0 \
    dash-table==4.6.2 \
    dash==1.11.0 \
    direct-redis \
    kaleido \
    matplotlib \
    numpy==1.19.5 \
    pandas>=1.3.0 \
    plotly \
    plotnine>=0.8.0 \
    psutil \
    pyarrow \
    seaborn==0.11.1 && \
    python3.8 -m pip install \
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
RUN python3.8 -m pip install \
    git+https://github.com/FredHutch/menu-driven-figure.git@e29bef4
ADD app/ /usr/bin/local/
RUN ln -s /usr/bin/local/gig-map /usr/local/bin/ && \
    ln -s /usr/bin/local/gig-map-cli /usr/local/bin/ && \
    gig-map --help && \
    gig-map-cli --help
RUN mkdir /work
WORKDIR /work