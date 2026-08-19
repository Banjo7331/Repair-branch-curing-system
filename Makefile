PY ?= python3
export PYTHONPATH := src

.DEFAULT_GOAL := help

.PHONY: help
help: ## Show this help
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

.PHONY: test
test: ## Run the test suite (mechanism, modality, annotation, reanalysis)
	$(PY) -m pytest tests

.PHONY: reference
reference: ## Re-run all three reference sets: mechanisms, modalities, episodes
	$(PY) -m repairbench.cli reference
	@echo
	$(PY) -m repairbench.cli reference --modalities
	@echo
	$(PY) -m repairbench.cli reference --reanalysis

.PHONY: annotation
annotation: ## Summarise the test annotation fixture
	$(PY) -m repairbench.cli annotation tests/data/mini.gff3

.PHONY: rules
rules: ## Print every rule, what it claims, and what it cites
	$(PY) -m repairbench.cli rules

.PHONY: lint
lint: ## Lint and type-check. Fails the build, like CI does.
	$(PY) -m ruff check src tests
	$(PY) -m mypy

.PHONY: check
check: lint test reference ## Everything CI runs
