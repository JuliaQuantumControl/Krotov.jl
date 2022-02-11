.PHONY: help test docs clean distclean devrepl codestyle
.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
    match = re.match(r'^([a-z0-9A-Z_-]+):.*?## (.*)$$', line)
    if match:
        target, help = match.groups()
        print("%-20s %s" % (target, help))
print("""
Instead of "make test", consider "make devrepl" if you want to run the test
suite or generate the docs repeatedly.

Make sure you have Revise.jl installed in your standard Julia environment
""")
endef
export PRINT_HELP_PYSCRIPT

GITORIGIN := $(shell git config remote.origin.url)

help:  ## show this help
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

JULIA ?= julia


test:  test/Manifest.toml ## Run the test suite
	$(JULIA) --project=test --color=auto --startup-file=yes --code-coverage="user" --depwarn="yes" --check-bounds="yes" -e 'include("test/runtests.jl")'
	@echo "Done. Consider using 'make devrepl'"


test/Manifest.toml: test/Project.toml ../scripts/installorg.jl
	$(JULIA) --project=test ../scripts/installorg.jl


docs/Manifest.toml: test/Manifest.toml
	cp test/*.toml docs/


devrepl: test/Manifest.toml docs/Manifest.toml ## Start an interactive REPL for testing and building documentation
	@$(JULIA) --project=test --banner=no --startup-file=yes -e 'include("test/init.jl")' -i


docs: docs/Manifest.toml  ## Build the documentation
	$(JULIA) --project=test docs/make.jl
	@echo "Done. Consider using 'make devrepl'"


clean: ## Clean up build/doc/testing artifacts
	rm -f src/*.cov test/*.cov examples/*.cov
	rm -f test/examples/*.*
	for file in examples/*.jl; do rm -f docs/src/"$${file%.jl}".*; done
	rm -rf docs/build


codestyle: test/Manifest.toml ../.JuliaFormatter.toml ## Apply the codestyle to the entire project
	$(JULIA) --project=test -e 'using JuliaFormatter; format(".", verbose=true)'
	@echo "Done. Consider using 'make devrepl'"


distclean: clean ## Restore to a clean checkout state
	rm -f Manifest.toml test/Manifest.toml
	rm -f docs/Manifest docs/Project.toml
	rm -rf test/data/* docs/data/*
	rm -rf test/examples/plots docs/src/examples/plots
