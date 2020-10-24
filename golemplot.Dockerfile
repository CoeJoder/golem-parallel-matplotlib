FROM python:3.6-slim

RUN apt update

RUN mkdir --parents /golem/work
COPY requirements.txt /golem/work/requirements.txt

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r /golem/work/requirements.txt

WORKDIR /golem/work
VOLUME /golem/work
