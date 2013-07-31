class GitError(RuntimeError):
    def __init__(self, exit_code):
        self.exit_code = exit_code
        RuntimeError.__init__(self, "git returned with non-zero exit code (%s)"%exit_code)

class DetachedHeadError(RuntimeError):
    def __init__(self):
        RuntimeError.__init__(self, "unexpectedly, git is in a detached HEAD state")
