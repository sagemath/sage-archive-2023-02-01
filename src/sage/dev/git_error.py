class GitError(RuntimeError):
    def __init__(self, exit_code):
        RuntimeError("git returned with non-zero exit code (%s)"%exit_code)

