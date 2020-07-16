from distutils import log
from distutils.command.build import build

class sage_build(build):
    sub_commands = [('build_cython', lambda *args: True)] + build.sub_commands

    def run_autogen(self):
        """
        Generate auto-generated sources.

        This must be done before building the python modules,
        see :trac:`22106`.
        """
        from sage_setup.autogen import autogen_all
        log.info("Generating auto-generated sources")
        for pkg in autogen_all():
            if pkg not in self.distribution.packages:
                    self.distribution.packages.append(pkg)

    def run(self):
        self.run_autogen()
        build.run(self)
