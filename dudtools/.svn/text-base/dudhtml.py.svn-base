"""Classes to generate HTML page with or without menu over the dud targets.

Michael Mysinger 200704 Created
"""

import os

STYLESHEET = "properties.css"

class BaseHTML:
    """Class to represent a simple html page without a menu."""
    def __init__(self, html, stylesheet, title):
        self.html = html
        self.stylesheet = stylesheet
        self.title = title

    def header(self):
        """Generate html header."""
        self.html.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n')
        self.html.write('<html>\n')
        self.html.write('<head>\n')
        self.html.write('    <title>%s</title>\n' % self.title)
        self.html.write('    <link rel="stylesheet" type="text/css" ' + 
                                 'href="%s" />\n' % self.stylesheet)
        self.html.write('</head>\n')
        self.html.write('<body>\n')
        self.html.write('<div id="whole_page">\n')

    def footer(self):
        """Generate html footer."""
        self.html.write('</div>\n')
        self.html.write('</body>\n')
        self.html.write('</html>\n')

    def content(self):
        """Default main html page contents."""
        self.html.write('<h1>%s</h1>\n' % self.title)
        self.html.write('<h2> Override content() to replace.</h2>\n')

    def main(self):
        """Generate main section."""
        self.html.write('<div id="main_page">\n')
        self.content()
        self.html.write('</div>\n')

    def menu(self):
        """Default BaseHTML is without a menu."""
        pass

    def page(self):
        """Generate full html page."""
        self.header()
        self.menu()
        self.main()
        self.footer()

class SingleThumbHTML(BaseHTML):
    """Show a thumbnail gallery of images on one page."""
    def __init__(self, html, stylesheet, title, images):
        BaseHTML.__init__(self, html, stylesheet, title)
        self.images = images

    def content(self):
        self.html.write('<h1>%s</h1>\n' % self.title)        
        self.html.write('<div id="image_menu">\n')
        for image in self.images:
            self.html.write('<a id="thumbnail" href="#nogo">\n') 
            self.html.write('    <img src="%s_thumbnail.png" />\n' % image)
            self.html.write('    <span><img src="%s.png" /></span>\n' % image)
            self.html.write('</a><br />\n')
        self.html.write('</div>\n')

class DudHTML(BaseHTML):
    """Class to represent a generic DUD html page with a menu."""
    def __init__(self, html, stylesheet, title, targets):
        BaseHTML.__init__(self, html, stylesheet, title)
        self.targets = targets

    def menu(self):
        """Generate menu of targets."""
        self.html.write('<div id="main_menu">\n')
        self.html.write('<div id="menu_item">\n')
        for target in self.targets:
            self.html.write('    <a href="%s.html">%s</a>\n' %
                       (target[1], target[0].upper()))
        self.html.write('</div>\n')
        self.html.write('</div>\n')

class ThumbHTML(DudHTML):
    """Show a thumbnail gallery of images for each target page."""
    def __init__(self, html, stylesheet, title, targets, images):
        DudHTML.__init__(self, html, stylesheet, title, targets)
        self.images = images

    def content(self):
        self.html.write('<h1>%s</h1>\n' % self.title)        
        self.html.write('<div id="image_menu">\n')
        for image in self.images:
            self.html.write('<a id="thumbnail" href="#nogo">\n') 
            self.html.write('    <img src="%s_thumbnail.png" />\n' % image)
            self.html.write('    <span><img src="%s.png" /></span>\n' % image)
            self.html.write('</a><br />\n')
        self.html.write('</div>\n')

class ImageHTML(DudHTML):
    """Quickly flip through one image per DUD target."""
    def __init__(self, html, stylesheet, title, targets, images):
        DudHTML.__init__(self, html, stylesheet, title, targets)
        self.images = images

    def content(self):
        self.html.write('<h1>%s</h1>\n' % self.title)        

    def menu(self):
        self.html.write('<div id="main_menu">\n')
        self.html.write('<div id="menu_item">\n')
        for target, image in zip(self.targets, self.images):
            self.html.write('<a id="quick" href="#nogo">%s\n'
                            % target[0].upper())
            self.html.write('    <img src="%s.png" />\n' % image)
            self.html.write('</a>\n')
        self.html.write('</div>\n')
        self.html.write('</div>\n')

