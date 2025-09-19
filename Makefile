.SHELLFLAGS  := -ilc
.ONESHELL:

ENV_NAME := nextflow
ENV_FILE := envs/nextflow.yaml
RUN_CONFIG ?= config/run01.config
RESOURCES_CONFIG := config/resources.config
RUN_ARGS ?= 
CLEAN_ARGS ?= -n
WHICH_FAILED = awk -F'\t' '!($$4 ~ /^[[:space:]]*OK[[:space:]]*$$/) { print $$3 }'

.PHONY: env_create env_remove nf_clean nf_run 

env_create:
	conda env create -n $(ENV_NAME) -f $(ENV_FILE)

env_remove:
	conda env remove -n $(ENV_NAME)

nf_clean:
	# Remove cached Nextflow files/dirs from failed runs (no OK in nextflow log)
	to_delete=$$(conda run -n $(ENV_NAME) nextflow log | $(WHICH_FAILED)); \
	test -z "$$to_delete" && { echo "Nothing to clean." >&2; exit 1; }; \
	echo "$$to_delete" | while IFS= read -r run_name; do \
		conda run -n $(ENV_NAME) nextflow clean "$$run_name" $(CLEAN_ARGS); \
	done

nf_run: 
	source ~/.bash_profile
	conda run -n $(ENV_NAME) \
		nextflow run . \
			-c $(RUN_CONFIG) \
			-c $(RESOURCES_CONFIG) \
			$(RUN_ARGS)
