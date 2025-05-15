
# exp3_rad_offset_gain

An example Python project focused on applying radiometric offset and gain corrections to satellite data. It includes utilities, a command-line interface, and testing & integration tools to work with satellite Level 1B data.

---

## What's Inside

### Structure

```aiignore
.
├── Dockerfile
├── Makefile
├── README.md
├── dev-requirements.in
├── dev-requirements.txt
├── integration
│   ├── nc_generator.py
│   ├── sample_l1b.nc
│   └── test_container.sh
├── mypy.ini
├── notebooks
│   ├── exp3.0_rad_offset.ipynb
│   └── settings
│       ├── gm_config.yaml
│       ├── l1bl2_config.yaml
│       ├── sgmgeo_config.yaml
│       ├── sgmrad_config.yaml
│       └── siml1b_config.yaml
├── pyproject.toml
├── rad_offset_gain*
├── requirements.in
├── requirements.txt
├── src
│   └── rad_offset_gain
│       ├── __init__.py
│       ├── cli.py
│       └── lib.py
├── tests
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_cli.py
│   └── test_lib.py
└── tox.ini

```

### Directory `src/rad_offset_gain/`

The core module of the project:
- `lib.py`: Implements the main logic for gain and offset correction.
- `cli.py`: Command-line interface for running the module.

### Directory `integration/`

Contains test data and scripts for validating the module in a near-real usage environment:
- `sample_l1b.nc`: Sample NetCDF files used for testing.
- `nc_generator.py`: Script to generate sample NetCDF files.
- `test_container.sh`: A shell script to test the Dockerized version of the project.

#### How to run?

```bash
$ sh integration/test_container.sh 
2025-05-11 14:00:48,924 DEBUG Adding radiance offset of 0.1 to sample_l1b.nc
2025-05-11 14:00:48,924 DEBUG Reading L1B data from sample_l1b.nc
2025-05-11 14:00:49,009 INFO Radiometric offset added successfully 
```

### Directory `tests/`

Unit tests written in `pytest` format to ensure the module works correctly:
- `conftest.py`: Pytest fixtures.
- `test_cli.py`: Tests for the command-line interface.
- `test_lib.py`: Tests for the core library functions.

#### How to run?

```bash
$ make tests
$ make coverage
```

### Directory `notebooks/`

Contains Jupyter notebooks that serve as exploratory tools:
- Demonstrate how the module works
- Provide visualization or exploratory testing

---

## Makefile Commands

The project is automation-friendly via a `Makefile`. You can use common development tasks with one command:

```sh
make venv                 # Create and prepare virtual environment
make install              # Install dependencies via tox
make uninstall            # Uninstall all project packages
make update-requirements  # Update requirements.txt from requirements.in
make fmt                  # Format code using Black
make lint                 # Lint code using Ruff
make mypy                 # Static type checking with MyPy
make flake8               # Lint using Flake8
make tests                # Run tests using tox (Python 3.10)
make coverage             # Run tests with coverage report
make build                # Build the project package
make docker               # Build a Docker image
make docker-bash          # Run Docker container in interactive bash
make container            # Build Docker image using Podman
make test-container       # Run integration tests in the container
make clean                # Remove temporary files and caches
make dist-clean           # Full cleanup, including virtual envs and caches
```

---

## Dependencies

Dependencies are managed via `tox` and defined in:
- `requirements.in`
  - compiled to `requirements.txt` using `pip-tools` (via `make update-requirements`)
- `dev-requirements.in`
  - compiled to `dev-requirements.txt` using `pip-tools` (via `make update-requirements`)
- `pyproject.toml`: Project configuration
- `mypy.ini`, `.flake8`: Linter and type-checking configs
- `tox.ini`: Tox configuration (managed via `Makefile`)

Main tools used:
- **Black** – code formatting
- **Ruff** – fast linter
- **MyPy** – static type checker
- **Tox** – testing and environment automation
- **Docker** and **Podman** – containerization

---

## Getting Started

### Running the CLI from source code

If you want to run the project locally directly from the source code, 
you can just run `rad_offset_gain` script from the command line.:

```bash
$ ./rad_offset_gain --help
usage: rad_offset_gain [-h] --input INPUT --output OUTPUT [--offset OFFSET]

Apply a radiance offset to L1B NetCDF data.

options:
  -h, --help       show this help message and exit
  --input INPUT    Input NetCDF file path
  --output OUTPUT  Output NetCDF file path
  --offset OFFSET  Offset fraction (e.g., 0.01 for 1%)
```

### Running the CLI from a system-wide installation

If you want to install the package system-wide, you can use `pip`:

```bash
# Test and build the package
$ make build

# Install the package
$ pip install dist/rad_offset_gain-1.0.0-py3-none-any.whl

# Run the CLI
$ rad_offset_gain --help
usage: rad_offset_gain [-h] --input INPUT --output OUTPUT [--offset OFFSET]
```

### Running the CLI from a virtual environment

If you want to install the package in a virtual environment:

```bash
# Create virtual environment, test and build the package
$ make build

# Run virtual environment and install the package
$ source .venv/bin/activate
(venv)$ pip install dist/rad_offset_gain-1.0.0-py3-none-any.whl

# Run the CLI
(venv)$ rad_offset_gain
usage: rad_offset_gain [-h] --input INPUT --output OUTPUT [--offset OFFSET]
```

### Running the CLI from a Docker container
If you want to run the package in a Docker container, you can use the provided `Dockerfile`:

```bash
# Build the Docker image
$ make docker

# Run the Docker container
$ docker run -it rad-offset-gain:latest --help
usage: cli.py [-h] --input INPUT --output OUTPUT [--offset OFFSET]

Apply a radiance offset to L1B NetCDF data.

options:
  -h, --help       show this help message and exit
  --input INPUT    Input NetCDF file path
  --output OUTPUT  Output NetCDF file path
  --offset OFFSET  Offset fraction (e.g., 0.01 for 1%)
```

See also Docker support section for more details.

---

## Running Tests

```bash
$ make tests
$ make coverage
```

---

## Docker Support

To test the project in a local containerized environment, you can use Docker. 
The `Dockerfile` is set up to build the project and run tests.

```bash
$ make docker         # Build image
```

In CI/CD pipelines, you can use the following commands to build and run the tests in a containerized environment. 
This command uses Podman instead of Docker.

```bash
$ make container      # Build image
$ make test-container # Run integration tests in the container
```

---

## Notebooks

Exploratory notebooks can be found in the `notebooks/` folder. These are useful for:
- Visualizing intermediate outputs
- Experimenting with different gain/offset values
- Manual testing and analysis

---

## Cleanup

```bash
make clean          # Remove temp files
make dist-clean     # Full project cleanup
```

---

## CI/CD

The workflow templates for these CI/CD pipelines were designed to be generic, so developers can easily test them locally using the `Makefile` - just like the pipeline would. This also makes it simple to adapt the templates to other modules in the project, improving scalability and portability.

The project uses GitHub Actions for CI/CD. 
There are two main workflows:
- **Python Package**: Runs on every push to the `main` branch and on pull requests. It includes:
  - Linting
  - Type checking
  - Unit tests
  - Coverage report
  - Pushing to PyPI _(commented out for now)_
- **Container**: Runs on every push to the `main` branch and on pull requests. It includes:
  - Building Python package
  - Building the Docker image
  - Running integration tests in the container
  - Pushing the Docker image to the registry _(commented out for now)_


### Examples

These examples utilize Self-Hosted runner deployed on GCP x64 Linux instance.
See details: https://github.com/undomained/teds/actions/runners?tab=self-hosted

#### Python Package

Example of a GitHub Actions workflow file for the Python package:

```yaml
# Template for building Python packages

name: Python Package Builder

on:
  workflow_call:
    inputs:
      project-path:
        required: true
        type: string

jobs:
  package-build:
    runs-on: self-hosted

    steps:
      - uses: actions/checkout@v4

      - name: Build release distributions
        run: |
          cd ${{ inputs.project-path }}
          make build

      - name: Upload distributions
        uses: actions/upload-artifact@v4
        with:
          name: release-dists
          path: ${{ inputs.project-path }}/dist/
```

#### Implementation of the Python package

- See [rad_offset_gain+python-publish.yml](.github/workflows/rad_offset_gain+python.yml) file for the Python package workflow source code.
- See [Actions/rad_offset_gain - Python Package](https://github.com/undomained/teds/actions/workflows/rad_offset_gain+python-publish.yml) for the latest run.

### Container

Example of a GitHub Actions workflow file for the container:

```yaml
# Template for building containers

name: Container Builder

on:
  workflow_call:
    inputs:
      project-path:
        required: true
        type: string

jobs:
  container-build:
    runs-on: self-hosted

    steps:
      - uses: actions/checkout@v4

      - name: Build Container
        run: |
          cd ${{ inputs.project-path }}
          make container

      - name: Test Container
        run: |
          cd ${{ inputs.project-path }}
          make test-container
```

#### Implementation of the container

- See [rad_offset_gain+python-publish.yml](.github/workflows/rad_offset_gain+container.yml) file for the container workflow source code.
- See [Actions/rad_offset_gain - Container Image](https://github.com/undomained/teds/actions/workflows/rad_offset_gain+container.yml) for the latest run.

### Notes

To keep costs low and shorten pipeline execution time, certain sections of the templates have been commented out, as this project serves only as a demonstration.

