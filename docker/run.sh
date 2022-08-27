#!/usr/bin/bash
if [[ "$#" -lt 1 ]]; then
  echo "Usage: $0 [args...]"
  exit 1
fi

name=$(python -c 'import uuid; print(str(uuid.uuid4()))')
echo "Container name: ${name}"

docker run -it \
  --name "${name}" \
  -v $(pwd):/host \
  --user "$(id -u):$(id -g)" \
  "vitreusx/pas-cg:latest" \
  $@

docker rm "${name}" 1>/dev/null