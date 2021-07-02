#sleep 86400
podman build -t zumis:v2.9.7 .
podman tag zumis:v2.9.7 tomkellygenetics/zumis:v2.9.7
podman tag zumis:v2.9.7 tomkellygenetics/zumis:latest
podman push tomkellygenetics/zumis:v2.9.7
podman push tomkellygenetics/zumis:latest
