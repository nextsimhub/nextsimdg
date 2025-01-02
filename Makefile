all: build

DOCKERFILE 	= Dockerfiles/Dockerfile.devenv

DOCKER_CMD	= docker
NAMESPACE 	= ghcr.io/nextsimhub
IMAGE_NAME	= nextsimdg-dev-env
TAG 		= latest
IMAGE		= $(NAMESPACE)/$(IMAGE_NAME):$(TAG)
OPTIONS	=
# OPTIONS	= --no-cache

.PHONY: build run

pull:
	$(DOCKER_CMD) pull $(IMAGE)

build:
	$(DOCKER_CMD) build -f $(DOCKERFILE) -t $(IMAGE) $(OPTIONS) .

run:
	$(DOCKER_CMD) run --rm -it -v /home/joe:/home/joe $(IMAGE)

push:
	$(DOCKER_CMD) push $(IMAGE)
