import fileinput






TAGS = [tag + '-TAG' for tag in ['FILTERED', 'DEREPLICATED', 'FIGURES', 'WARNINGS']]

class HTMLBuilder():

    def __init__(self, html_path):
        self.html_path = html_path
        self.replacements = dict()


    def build_filtered(self):

        pass


    def build_pdfs(self, pdfs):

        def _pdfTag(pdf):
            return f'<embed src="figures/{pdf}" width="1000px" height="1000px">'

        rep = ''

        for pdf in pdfs:
            rep += self._encase_p(pdf) + self._encase_p(_pdfTag(pdf)) + '\n'

        self.replacements['FIGURES-TAG'] = rep


    def build_warnings(self, warnings):
        if warnings.strip() == '':
            warnings = 'No warnings about almost divergent secondary clustering or remaining genome similarity'

        self.replacements['WARNINGS-TAG'] = warnings


    def build(self):

        
        for line in fileinput.input(self.html_path, inplace=True):
            line_stripped = line.strip()
            if line_stripped in TAGS and line_stripped in self.replacements:
                print(self.replacements[line_stripped])
            else:
                print(line, end='')

    def _encase_p(self, paragraph):
        return '<p>' + paragraph + '</p>'



















