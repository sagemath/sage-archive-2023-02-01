ARG BASE_GITHUB_REPOSITORY=sagemath/sage
ARG BASE_TAG=dev
FROM ghcr.io/${BASE_GITHUB_REPOSITORY}/sage-docker-gitpod-standard-with-targets:${BASE_TAG} as with-targets
# Fast doc rebuilds do not work because
# "loading pickled environment... failed; source directory has changed"
# So remove the intermediate output to save space (1GB), but keep the built HTML documentation.
RUN (cd /home/gitpod/sage-local/share/doc/sage && rm -Rf doctrees inventory)

FROM ghcr.io/${BASE_GITHUB_REPOSITORY}/sage-docker-gitpod-standard-with-system-packages:${BASE_TAG}
COPY --chown=gitpod:gitpod --from=with-targets /home/gitpod/sage/logs  /home/gitpod/sage/logs
COPY --chown=gitpod:gitpod --from=with-targets /home/gitpod/sage-local /home/gitpod/sage-local
