"""The iCH360 site provides drawings without SBML model downloads."""
from pathlib import Path
import tempfile
import unittest
from prepare_github_pages import check_ich360_publication


class PublicationTests(unittest.TestCase):
    def test_published_ich360_contains_no_model_files_or_links(self):
        check_ich360_publication(Path(__file__).resolve().parent / 'docs')

    def test_rejects_source_export_and_compressed_model_files(self):
        for name in ('model.sbml', 'source.xml', 'layout.SBML', 'model.sbml.gz'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp) / 'eco_ich_360'
                folder.mkdir()
                (folder / name).write_text('<sbml/>')
                with self.assertRaisesRegex(ValueError, 'must not be published'):
                    check_ich360_publication(Path(tmp))

    def test_rejects_local_external_and_encoded_model_links(self):
        for link in ('layout.sbml', 'https://example.org/ich.xml?download=1', 'layout%2Esbml'):
            with self.subTest(link=link), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp) / 'eco_ich_360'
                folder.mkdir()
                (folder / 'index.html').write_text(f'<a href="{link}">Model</a>')
                with self.assertRaisesRegex(ValueError, 'must not link'):
                    check_ich360_publication(Path(tmp))
