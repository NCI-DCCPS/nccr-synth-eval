"""Boundary and CLI checks for the standalone three-CSV evaluator."""
from collections import Counter
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'compare_iccc_age_csv.py'
spec = importlib.util.spec_from_file_location('csv_eval', SCRIPT)
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)


class CsvEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def csv(self, name, rows, header=('icccMajor','ageRecodeSingle90p')):
        path = self.root/name
        with path.open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        return path

    def run_cli(self, paths, extra=(), html=False):
        args = [str(p) for p in paths]+['--out',str(self.root/'report')]+list(extra)
        if not html:
            args.append('--no-html')
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            e.main(args)
        return json.loads((self.root/'report/evaluation.json').read_text(encoding='utf-8'))

    def test_identical_inputs_normalization_blanks_and_unknown(self):
        p = self.csv('data.csv',[('1','0'),('02','99'),('',''),('12','90')])
        r = self.run_cli([p,p,p])
        self.assertEqual(r['inputs']['Reference']['included_rows'],4)
        self.assertEqual(r['inputs']['Reference']['code_normalization_counts']['icccMajor'],1)
        self.assertTrue(all(m['joint_tv']==0 for m in r['metrics']))
        text = (self.root/'report/cell_comparison.csv').read_text(encoding='utf-8')
        self.assertIn('Blank / missing',text)
        self.assertIn('99',text)

    def test_filter_custom_columns_and_empty_index_subset(self):
        p = self.csv('data.csv',[('1','4','2000','1'),('2','5','2023','1'),('1','3','','1')],
                     ('major','age','year','index'))
        r = self.run_cli([p,p,p],['--iccc-column','major','--age-column','age','--year-column','year',
                                 '--index-column','index','--min-year','1995','--max-year','2022'])
        self.assertEqual(r['inputs']['Reference']['included_rows'],1)
        self.assertEqual(r['inputs']['Reference']['excluded_year_rows'],2)
        self.assertEqual(r['inputs']['Reference']['invalid_or_missing_year_rows'],1)
        self.assertEqual(len(r['skipped_subsets']),4)

    def test_no_shared_support_html(self):
        a = self.csv('a.csv',[('01','00')])
        b = self.csv('b.csv',[('02','10')])
        r = self.run_cli([a,b,b],html=True)
        self.assertTrue(all(m['iccc_given_age_weighted_tv'] is None for m in r['metrics']))
        text = (self.root/'report/report.html').read_text(encoding='utf-8')
        self.assertIn('N/A',text)
        self.assertNotIn('frozen G1 reference',text)
        self.assertEqual(text.count('class="plotly-graph-div"'),4)

    def test_invalid_age_and_group_disagreement(self):
        for age, band in [('91','19'),('05','01')]:
            with self.subTest(age=age):
                p = self.csv('a.csv',[('01',age,band)],('icccMajor','ageRecodeSingle90p','ageRecode20Groups'))
                with self.assertRaises(ValueError):
                    self.run_cli([p,p,p])
                self.assertFalse((self.root/'report').exists())

    def test_missing_and_duplicate_headers(self):
        for header in [('icccMajor','other'),('icccMajor','icccMajor')]:
            with self.subTest(header=header):
                p = self.csv('a.csv',[('01','00')],header)
                with self.assertRaises(ValueError):
                    self.run_cli([p,p,p])

    def test_existing_output_refused_and_cli_failure_is_nonzero(self):
        p = self.csv('a.csv',[('01','00')])
        out = self.root/'report'
        out.mkdir()
        sentinel = out/'sentinel'
        sentinel.write_text('preserve')
        proc = subprocess.run([sys.executable,str(SCRIPT),str(p),str(p),str(p),'--out',str(out),'--no-html'],
                              capture_output=True,text=True)
        self.assertNotEqual(proc.returncode,0)
        self.assertEqual(sentinel.read_text(),'preserve')

    def test_conditional_tv_separates_age_marginal_and_support(self):
        r = Counter({('01','00'):90,('01','01'):10})
        g = Counter({('01','00'):10,('01','01'):90})
        m = e.summarize_counts(r,g)[0]
        self.assertAlmostEqual(m['joint_tv'],.8)
        self.assertEqual(m['iccc_given_age_weighted_tv'],0)
        m = e.summarize_counts(Counter({('01','00'):100}),Counter({('01','00'):50,('02','99'):50}))[0]
        self.assertEqual(m['conditional_coverage_fraction'],.5)
        self.assertEqual(m['generated_rows_in_reference_absent_ages'],50)

    def test_unicode_labels_and_version(self):
        p = self.csv('data.csv',[('01','00')])
        r = self.run_cli([p,p,p],['--labels','生成 A','生成 B'],html=True)
        self.assertEqual(r['runtime']['evaluator'],'0.1.0')
        self.assertIn('生成 A',(self.root/'report/summary.md').read_text(encoding='utf-8'))
        self.assertIn('生成 B',(self.root/'report/report.html').read_text(encoding='utf-8'))
        proc = subprocess.run([sys.executable,str(SCRIPT),'--version'],capture_output=True,text=True)
        self.assertEqual(proc.returncode,0)
        self.assertEqual(proc.stdout.strip(),'iccc-age 0.1.0')


if __name__ == '__main__':
    unittest.main()
