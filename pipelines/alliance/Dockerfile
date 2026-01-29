FROM python:3.8-bullseye

WORKDIR /usr/src/gene_descriptions

# Copy and install dependencies first (better layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Install the genedescriptions package
RUN pip install --no-cache-dir -e .

# Create output directory
RUN mkdir -p pipelines/alliance/generated_descriptions

# Run the alliance pipeline in parallel mode
CMD ["python", "-m", "pipelines.alliance.alliance_pipeline", \
     "-c", "pipelines/alliance/config_alliance.yml", \
     "-L", "INFO", \
     "--parallel"]