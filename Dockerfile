FROM python:3.8-slim as build

WORKDIR /usr/src/app

ENV PATH=/root/.local/bin:$PATH PYTHONUNBUFFERED=1

RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential gcc
RUN pip install --upgrade pip wheel

COPY setup.py .
COPY requirements.txt .
COPY README.md .
RUN python setup.py deps

COPY . .
RUN pip install --user .

FROM python:3.8-slim as main

ENV PATH=/root/.local/bin:$PATH PYTHONUNBUFFERED=1
ENV PYTHONPATH=/root/.local/lib/python3.8/site-packages:$PYTHONPATH
ENV LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib

COPY --from=build /root/.local /root/.local

CMD ['simulate-cancer-cells']
CMD ['optimize-vaccine-ilp']
CMD ['create-bar-chart']
CMD ['evaluate-vaccine-response']
