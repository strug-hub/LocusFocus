FROM python:3.10-slim as base

ARG USERNAME=flask
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN mkdir code \
  && groupadd --gid $USER_GID $USERNAME \
  && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
  && chown -R ${USER_UID}:${USER_GID} /code

RUN apt-get update && \
  apt-get install -y \
  build-essential \
  curl \
  gfortran \
  libblas-dev \
  liblapack-dev \
  libpcre2-dev \
  r-base \
  r-base-dev \
  xauth \
  vim \
  wget

# install plink
RUN wget -O /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
  && unzip -d /tmp /tmp/plink.zip \
  && mv /tmp/plink /usr/local/bin/plink \
  && rm /tmp/*

# install liftover
RUN curl -f -L -O https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver \
  -f -L -O https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz \
  -f -L -O https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz && \
  chmod +x ./liftOver && \
  mv liftOver /usr/local/bin/liftOver && \
  mkdir /usr/local/share/liftOver && \
  mv hg38ToHg19.over.chain.gz hg19ToHg38.over.chain.gz /usr/local/share/liftOver/

# Install Poetry
# https://github.com/python-poetry/poetry/issues/6397#issuecomment-1236327500
ENV POETRY_HOME=/opt/poetry

RUN python3 -m venv $POETRY_HOME

RUN $POETRY_HOME/bin/pip install -U pip \
  && $POETRY_HOME/bin/pip install poetry==1.5.1

ENV VIRTUAL_ENV=/poetry-env \
  PATH="/poetry-env/bin:$POETRY_HOME/bin:$PATH"

RUN python3 -m venv $VIRTUAL_ENV \
  && chown -R $USER_UID:$USER_GID $POETRY_HOME /poetry-env

WORKDIR /code

# handle R libraries

ENV R_LIBS_USER=/home/${USERNAME}/Rlibs

# only root can write here, so we'll save these into the default package dir and use to store application packages in the user's dir
RUN Rscript -e "install.packages(c('littler', 'docopt'))" \
  && ln -s /usr/local/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
  && ln -s /usr/local/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
  && ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/r

USER $USERNAME

RUN mkdir /home/${USERNAME}/Rlibs \
  && echo "R_LIBS=/home/${USERNAME}/Rlibs" > /home/${USERNAME}/.Renviron

# Install R packages (order matters, evidently)
RUN install2.r --error \
  stringr \
  data.table \
  argparser \
  CompQuadForm \
  here \
  BiocManager \
  remotes \
  zeallot

RUN R -e "BiocManager::install('GenomicRanges')"
RUN R -e "BiocManager::install('biomaRt')"

# this seems not to respect user install path, might need to install as root and copy over
# RUN R -e "remotes::install_version('Matrix', version = '1.2')" # this seems to be already installed via other deps, however

# Link plink to work dir
RUN ln -s /usr/local/bin/plink /code/plink
RUN ln -s /usr/local/bin/liftOver /code/liftOver

FROM base AS dev

COPY --chown=$USERNAME:$USERNAME ./pyproject.toml /code/pyproject.toml
COPY --chown=$USERNAME:$USERNAME ./poetry.lock /code/poetry.lock
COPY --chown=$USERNAME:$USERNAME ./README.md /code/README.md
COPY --chown=$USERNAME:$USERNAME ./app /code/app
COPY --chown=$USERNAME:$USERNAME ./tests /code/app

RUN poetry install --with dev

FROM base AS prod

COPY --chown=$USERNAME:$USERNAME ./pyproject.toml /code/pyproject.toml
COPY --chown=$USERNAME:$USERNAME ./poetry.lock /code/poetry.lock
COPY --chown=$USERNAME:$USERNAME ./README.md /code/README.md
COPY --chown=$USERNAME:$USERNAME ./app /code/app

RUN poetry install --no-dev
