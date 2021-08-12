[build-system]
# Minimum requirements for the build system to execute.
requires = [
    'sage-conf',
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        cython         \
                    ')]
build-backend = "setuptools.build_meta"
