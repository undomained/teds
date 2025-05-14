# Use a Python base image
FROM python:3.10-slim

# Set environment variables
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install system dependencies
RUN apt-get update

# && apt-get install -y gcc libffi-dev && rm -rf /var/lib/apt/lists/*

# Create virtualenv and install tox
RUN python -m venv $VIRTUAL_ENV
RUN pip install --upgrade pip tox

# Set workdir
WORKDIR /app

# Copy project
COPY dist/ /app/dist/

RUN pip install /app/dist/rad_offset_gain-*.whl && rm -rf /app/dist

ENTRYPOINT ["python", "-m", "rad_offset_gain.cli"]

