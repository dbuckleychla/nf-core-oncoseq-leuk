"""Caller display filters preserve native PASS versus filtered evidence."""
import json
import shutil
import subprocess
import unittest
from pathlib import Path


@unittest.skipUnless(shutil.which('node'), 'Node is needed for report JavaScript checks')
class CallerFilterTests(unittest.TestCase):
    def test_pass_any_and_missing_evidence(self):
        template = (Path(__file__).resolve().parents[2] / 'assets/review/report.html').read_text()
        start = template.index('function matchesCallers(')
        end = template.index('function table(', start)
        script = template[start:end] + '''
        const row={evidence:[{caller:'clairsto',filter:['PASS']},{caller:'deepsomatic',filter:['LowQual']}]};
        const choices=(mode, callers)=>['clair3','clairsto','deepsomatic'].map(c=>({caller:c,mode,checked:callers.includes(c)}));
        const result={
          filteredPass:matchesCallers(row,choices('PASS',['deepsomatic']),'ANY'),
          filteredAny:matchesCallers(row,choices('ANY',['deepsomatic']),'ANY'),
          bothPass:matchesCallers(row,choices('PASS',['clairsto','deepsomatic']),'ALL'),
          eitherPass:matchesCallers(row,choices('PASS',['clairsto','deepsomatic']),'ANY'),
          unchecked:matchesCallers(row,choices('ANY',[]),'ANY'),
          defaultMissing:matchesCallers({evidence:[]},choices('ANY',['clair3','clairsto','deepsomatic']),'ANY'),
          mixedFilter:matchesCallers({evidence:[{caller:'clairsto',filter:['PASS','LowQual']}]},choices('PASS',['clairsto']),'ANY'),
          stringPass:matchesCallers({evidence:[{caller:'clairsto',filter:'PASS'}]},choices('PASS',['clairsto']),'ANY'),
          recordedAgreement:matchesCallers({evidence:[],exact_pass_agreement:true},choices('PASS',['clairsto','deepsomatic']),'ALL'),
        };
        process.stdout.write(JSON.stringify(result));
        '''
        result = subprocess.run(['node', '--input-type=module', '-e', script], check=True, capture_output=True, text=True)
        self.assertEqual(json.loads(result.stdout), dict(filteredPass=False, filteredAny=True, bothPass=False,
                         eitherPass=True, unchecked=False, defaultMissing=True, mixedFilter=False,
                         stringPass=True, recordedAgreement=True))

    def test_native_read_fraction_and_significance_categories(self):
        template = (Path(__file__).resolve().parents[2] / 'assets/review/report.html').read_text()
        start = template.index('function significanceCategory(')
        end = template.index('function cosmicCell(', start)
        script = "function annotationText(v){return Array.isArray(v)?v.join(', '):String(v||'');}\n" + template[start:end] + """
        console.log(JSON.stringify({complete:readSupport({supporting_reads:3,depth:20}),missing:readSupport({af:.5,depth:20}),
          empty:readSupport(),zero:readSupport({supporting_reads:0,depth:20}),
          conflict:significanceCategory('Conflicting_classifications_of_pathogenicity'),
          vus:significanceCategory('Uncertain_significance'),pathogenic:significanceCategory('Pathogenic'),benign:significanceCategory('Benign')}));
        """
        result=subprocess.run(['node','--input-type=module','-e',script],check=True,capture_output=True,text=True)
        self.assertEqual(json.loads(result.stdout),dict(complete='3/20',missing='—/20',empty='—/—',zero='0/20',
                         conflict='uncertain',vus='uncertain',pathogenic='pathogenic',benign='neutral'))
