#!/bin/bash

sudo docker image rm -f registry.hzdr.de/ufz/tb3-cite/cbc/mag2ena:latest
sudo docker run -p 8501:8501 registry.hzdr.de/ufz/tb3-cite/cbc/mag2ena:latest