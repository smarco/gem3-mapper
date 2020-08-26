# Build the image:
docker build -t gem-image ./ 

# Additional arguments: (pass them to the build command)
  GEM_MAPPER_VERSION
  INSTALL_BASE
  SRC_BASE

docker build -t gem-image ./ --build-arg GEM_MAPPER_VERSION=v3.3.0

