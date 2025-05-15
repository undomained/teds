# Makefile for rad_offset_gain module

PYTHON = .venv/bin/python
PIP = .venv/bin/pip
TOX = .venv/bin/tox
BLACK = .venv/bin/black
RUFF = .venv/bin/ruff

.PHONY: update-requirements fmt install uninstall build tests coverage clean flake8

all: help

dist: all

help:
	@echo "Makefile for rad_offset_gain module"
	@echo "Available targets:"
	@echo "  all                 - Show this help message"
	@echo "  build               - Build distribution packages"
	@echo "  clean               - Clean build artifacts and cache"
	@echo "  coverage            - Run tests with coverage"
	@echo "  dist                - Run all targets and build distribution packages"
	@echo "  dist-clean          - Clean all build artifacts and caches"
	@echo "  docker              - Build Docker image"
	@echo "  docker-bash         - Run Docker container with bash"
	@echo "  fmt                 - Format code using black"
	@echo "  install             - Install project dependencies"
	@echo "  lint                - Linting with ruff"
	@echo "  mypy                - Run MyPy type checking"
	@echo "  tests               - Run tests using tox"
	@echo "  uninstall           - Uninstall project and its dependencies"
	@echo "  update-requirements - Update requirements.txt from requirements.in"
	@echo "  venv                - Create a virtual environment and install dependencies"
	@echo ""
	@echo "  help                - Show this help message"
	@echo ""


# Virtual environment setup
.PHONY: venv
venv:
	bash -c "python3 -m venv .venv && \
	source .venv/bin/activate && \
	pip install --upgrade pip && \
	pip install tox"

# Update requirements.txt from requirements.in using pip-tools
update-requirements: venv
	$(TOX) -e update-requirements

# Format code using black
fmt: venv
	$(TOX) -e fmt

# Install project dependencies
install: venv
	$(TOX) -e install

# Uninstall project and its dependencies (except pip-tools and black)
uninstall:
	$(TOX) -e uninstall

# Build distribution packages
build: venv
	$(TOX) -e build

# Clean build artifacts and cache
clean:
	rm -rf build dist src/*.egg-info
	find . -type d -name '__pycache__' -exec rm -rf {} +
	find . -type f -name '*.pyc' -delete

dist-clean: clean
	rm -rf .tox .venv .pytest_cache .mypy_cache .coverage
	rm -rf .eggs .hypothesis .ruff_cache
	rm -rf coverage.xml htmlcov

# Run tests using tox
tests: venv
	$(TOX) -e py310 -- tests

# Run tests with coverage
coverage: venv
	$(TOX) -e py310 -- --cov=rad_offset_gain --cov-report=term-missing tests

# MyPy type checking
mypy: venv
	$(TOX) -e mypy

# Linting with ruff
lint: venv
	$(TOX) -e ruff

# Linting with flake8
flake8: venv
	$(TOX) -e flake8

# Build Docker image
docker: build
	$(TOX) -e docker

docker-bash:
	docker run -ti rad-offset-gain bash

# Build Docker image using podman
container: build
	$(TOX) -e container

test-container:
	$(TOX) -e test-container