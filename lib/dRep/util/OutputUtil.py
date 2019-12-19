import fileinput






TAGS = [tag + '-TAG' for tag in ['FILTERED', 'DEREPLICATED', 'FIGURES', 'WARNINGS']]

class HTMLBuilder():

    def __init__(self, html_path):
        self.html_path = html_path
        self.replacements = dict()


    def build_filtered(self):

        pass


    def build_pdfs(self, pdfs):

        def build_htmlPdfTag(pdf):
            return fr'<a href="img/{pdf}>{pdf}<\a>'

        rep = ''

        for pdf in pdfs:
            rep += build_htmlPdfTag(pdf) + '\n'

        self.replacements['FIGURES-TAG'] = rep


    def build_warnings(self, warnings):
        if warnings.strip() == '':
            warnings = 'No warnings'

        self.replacements['WARNINGS-TAG'] = warnings


    def build(self):

        
        for line in fileinput.input(self.html_path, inplace=True):
            for TAG in TAGS:
                if line.strip() == TAG:
                    print(self.replacements[TAG])





















