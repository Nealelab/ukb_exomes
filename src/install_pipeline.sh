#!/usr/bin/env bash


ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ "${ROLE}" == 'Master' ]]; then

    /opt/conda/bin/pip install --upgrade cerberus aiohttp jwt
    cd /home/hail
    git clone https://github.com/hail-is/hail.git
    git remote add jigold https://github.com/jigold/hail.git
    git fetch --all

    apt-get update
    apt-get -y install \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg2 \
        software-properties-common \
        tabix
    curl -fsSL https://download.docker.com/linux/debian/gpg | sudo apt-key add -
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
    apt-get update
    apt-get install -y --allow-unauthenticated docker-ce
    usermod -a -G docker konradk
fi
