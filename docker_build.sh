#sleep 86400
podman build -t zumis2:v2.9.7 .
podman tag zumis2:v2.9.7 chrzie/zumis2:v2.9.7
podman tag zumis2:v2.9.7 chrzie/zumis2:latest
podman push chrzie/zumis2:v2.9.7
podman push chrzie/zumis2:latest
