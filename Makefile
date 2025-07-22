.SHELLFLAGS  := -ilc
.ONESHELL:

ENV_NAME := nextflow
ENV_FILE := envs/nextflow.yaml
RUN_CONFIG ?= config/run01.config
RESOURCES_CONFIG := config/resources.config
ARGS ?= 


.PHONY: env clean run

env:
	conda env create -n $(ENV_NAME) -f $(ENV_FILE)

clean:
	conda env remove -n $(ENV_NAME)

run:            
	source ~/.bashrc
	source ~/.bash_profile
	conda run -n $(ENV_NAME) \
		nextflow run . \
			-c $(RUN_CONFIG) \
			-c $(RESOURCES_CONFIG) \
			$(ARGS)