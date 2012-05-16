# -*- coding: utf-8 -*-
"""
    inventory builder
    ~~~~~~~~~~~~~~~~~

    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
"""
from sphinx.builders.html import StandaloneHTMLBuilder
from sphinx.util.console import bold

class InventoryBuilder(StandaloneHTMLBuilder):
    """
    A customized HTML builder which only generates intersphinx "object.inv"
    inventory files and pickle files. The documentation files are not written.
    """
    name = "inventory"

    def write_doc(self, docname, doctree):
        """
        Don't write any doc
        """

    def finish(self):
        """
        Only write the inventory files.
        """
        self.write_buildinfo()
        self.dump_inventory()

    def removed_method_error(self):
        """
        Raise an error if this method is called.

        This is just for making sure that some writer methods are indeed
        deactivated.
        """
        raise RuntimeError, "This function shouldn't be called in \"%s\" builder"%(self.name)

    copy_image_files = removed_method_error
    copy_download_files = removed_method_error
    copy_static_files = removed_method_error
    handle_finish = removed_method_error

def setup(app):
    app.add_builder(InventoryBuilder)

