#!/bin/zsh
set -e
cd -- "${0:A:h}"
exec python3 -m spatial_lab.server --port 8841
