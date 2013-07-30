class GitError(RuntimeError):
    def __init__(self, exit_code):
        RuntimeError.__init__(self, "git returned with non-zero exit code (%s)"%exit_code)

