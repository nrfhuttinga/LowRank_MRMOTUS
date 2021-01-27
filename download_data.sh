#!/bin/bash

URL="https://surfdrive.surf.nl/files/index.php/s/QdOryuR8GKVSzao/download"

curl -o temp.zip $URL && unzip -o temp.zip && rm -f temp.zip
