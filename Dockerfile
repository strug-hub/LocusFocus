FROM python:3.10-buster as base

ARG USERNAME=flask
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN mkdir code \
  && groupadd --gid $USER_GID $USERNAME \
  && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
  && chown -R ${USER_UID}:${USER_GID} /code

# can upgrade R to something sensible later, once we have the package lists and versions set
# https://github.com/rocker-org/rocker/blob/2f92c6c8b8da7b3e61aabc44cacc0439bf267d31/r-base/3.6.3/Dockerfile

COPY --from=r-base:3.6.3 /usr/lib/R /usr/lib/R
COPY --from=r-base:3.6.3 /usr/local/lib/R /usr/local/lib/R
COPY --from=r-base:3.6.3 /usr/lib/libR.so /usr/lib/libR.so
COPY --from=r-base:3.6.3 /usr/bin/R /usr/bin/R
COPY --from=r-base:3.6.3 /usr/bin/r /usr/bin/r
COPY --from=r-base:3.6.3 /usr/bin/Rscript /usr/bin/Rscript
COPY --from=r-base:3.6.3 /lib/x86_64-linux-gnu/libreadline.so.8 /lib/x86_64-linux-gnu/libreadline.so.8
COPY --from=r-base:3.6.3 /usr/lib/x86_64-linux-gnu/libm.so /usr/lib/x86_64-linux-gnu/libm.so
COPY --from=r-base:3.6.3 /lib/x86_64-linux-gnu/libm.so.6 /lib/x86_64-linux-gnu/libm.so.6
COPY --from=r-base:3.6.3 /etc/R /etc/R
COPY --from=r-base:3.6.3 /usr/share/R /usr/share/R

RUN chown -R root:${USERNAME} /usr/local/lib/R

RUN ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
  && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r

ENV EDITOR vim
ENV R_BASE_VERSION 3.6.3
ENV R_LIBS_USER /home/${USERNAME}/Rlibs

RUN apt-get update && \
  apt-get install -y \
  libblas-dev \
  build-essential \
  gfortran \
  libblas-dev \
  liblapack-dev \
  libpcre2-dev \
  xauth \
  vim

# install plink
RUN wget -O /tmp/plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip \
  && unzip -d /tmp /tmp/plink.zip \
  && mv /tmp/plink /usr/local/bin/plink \
  && rm /tmp/*

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

USER $USERNAME

RUN mkdir /home/${USERNAME}/Rlibs \
  && echo "R_LIBS=/home/${USERNAME}/Rlibs" > /home/${USERNAME}/.Renviron

RUN R -e 'install.packages("docopt")'

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
