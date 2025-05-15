# Makefile for rad_offset_gain module

PYTHON = .venv/bin/python
PIP = .venv/bin/pip
TOX = .venv/bin/tox
BLACK = .venv/bin/black
RUFF = .venv/bin/ruff

.PHONY: update-requirements fmt install uninstall build
.PHONY: tests coverage clean flake8 mypy lint
.PHONY: dist-clean docker container test-container
.PHONY: all help venv

all: help

help:
	@echo "Makefile for rad_offset_gain module"
	@echo "Available targets:"
	@echo "  all                 - Show this help message"
	@echo "  build               - Build distribution packages"
	@echo "  clean               - Clean build artifacts and cache"
	@echo "  coverage            - Run tests with coverage"
	@echo "  dist-clean          - Clean all build artifacts and caches"
	@echo "  docker              - Build Docker image"
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

# Run all tests using tox
tests: venv
	$(TOX) -e run-tests

# Run tests with coverage
coverage: venv
	$(TOX) -e coverage

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

# Build Docker image using podman
container: build
	$(TOX) -e container

test-container:
	$(TOX) -e test-container

test-all-targets: clean dist-clean venv dist update-requirements fmt tests \
	coverage mypy lint flake8 build install uninstall docker dist-clean
