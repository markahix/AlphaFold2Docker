# AlphaFold2Docker

This is a very simple setup.  It's ColabFold, topped off with a Flask WebUI, and stuffed into a Docker container.

If you already have Docker or DockerDesktop set up on your system, you simply need to pull the image from [DockerHub](https://hub.docker.com/repository/docker/markahix/alphafold2).

    docker run --gpus all -p 5050:8543 markahix/alphafold2:latest

