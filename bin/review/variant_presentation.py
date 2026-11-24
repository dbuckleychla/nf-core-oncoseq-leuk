"""Explain existing selection and assign one display group per allele."""
import re


def text(value):
    if isinstance(value, (list, tuple)):
        value = ', '.join(str(v) for v in value if v is not None)
    return str(value or '').replace('_', ' ').replace('&', ', ').strip()


def significance(value, term):
    return bool(re.search(r'(^|[|/,])\s*' + term + r'\s*($|[|/,])', text(value).lower()))


def cosmic_summary(row):
    sources = [s for e in row.get('evidence', []) for s in e.get('cosmic', {}).get('sources', [])]
    identifiers = sorted({i for s in sources for i in s.get('identifiers', [])})
    targeted = {t['identifier']: t for e in row.get('evidence', [])
                for t in e.get('cosmic', {}).get('targeted_evidence', [])}
    samples, blood_samples, publications, statuses = set(), set(), set(), set()
    complete = bool(targeted) and set(targeted) == set(identifiers) and all(t.get('status') == 'matched' and 'evidence' in t for t in targeted.values())
    for target in targeted.values():
        blood_ids = {p['phenotype_id'] for p in target.get('phenotypes', [])
                     if (p.get('classification') or {}).get('PRIMARY_SITE') == 'haematopoietic_and_lymphoid_tissue'}
        for evidence in target.get('evidence', []):
            sample = evidence.get('COSMIC_SAMPLE_ID')
            if sample and sample not in ('.', 'NS'):
                samples.add(sample)
                if evidence.get('COSMIC_PHENOTYPE_ID') in blood_ids:
                    blood_samples.add(sample)
            pmid = evidence.get('PUBMED_PMID')
            if pmid and pmid not in ('.', 'NS'):
                publications.add(pmid)
            status = evidence.get('MUTATION_SOMATIC_STATUS')
            if status:
                statuses.add(status.lower())
    categories = []
    if any('confirmed somatic' in s for s in statuses):categories.append('Confirmed somatic reports')
    if any('another cancer sample' in s for s in statuses):categories.append('Other somatic reports')
    if any('unknown origin' in s for s in statuses):categories.append('Unknown origin reports')
    if statuses and not categories:categories.append('Other reported status')
    return dict(matched=bool(identifiers), identifiers=identifiers,
                tissue_known=any(p.get('classification') for t in targeted.values() for p in t.get('phenotypes', [])),
                blood_match=any(e.get('cosmic', {}).get('blood_match') for e in row.get('evidence', [])),
                sample_count=len(samples) if complete else None,
                blood_sample_count=len(blood_samples) if complete else None,
                publication_count=len(publications) if complete else None,
                somatic_categories=categories, counts_complete=complete)


def enrich_variants(data):
    for row in data.get('small_variants', {}).get('records', []):
        selections = row.get('selections', [])
        primary = min(selections, key=lambda s: (s.get('priority', 99), s.get('caller') == 'consensus'), default={})
        native_reasons = set(r for s in selections for r in [s.get('reason'), *s.get('inclusion_reasons', [])] if r)
        assertions = []
        for item in [*selections, *row.get('evidence', [])]:
            for assertion in item.get('clinvar', []):
                if assertion not in assertions:
                    assertions.append(assertion)
        sig = '|'.join(text(a.get('CLNSIG')) for a in assertions)
        conflict = 'conflict' in sig.lower() or any(a.get('CLNSIGCONF') or 'conflicting' in text(a.get('CLNREVSTAT')).lower() for a in assertions)
        benign = significance(sig, r'(?:likely )?benign')
        uncertain = significance(sig, r'uncertain significance')
        plp = significance(sig, r'(?:likely )?pathogenic')
        reviewed = any(significance(a.get('CLNSIG'), r'(?:likely )?pathogenic') and
                       any(k in text(a.get('CLNREVSTAT')).lower() for k in ('criteria provided', 'expert panel', 'practice guideline')) for a in assertions)
        clean = reviewed and not (conflict or benign or uncertain)
        lof = 'Predicted loss of function' in native_reasons and not (conflict or benign)
        panel_coding = any('Panel coding variant' in r for r in native_reasons)
        transcripts = []
        for selection in selections:
            candidates = [(selection.get('gene_context', {}).get('transcript') or {}, 'Selected COSMIC transcript' if selection.get('cosmic_policy') else 'Selected transcript')]
            candidates += [(t, 'VEP annotation' if selection.get('cosmic_policy') else 'Loss-of-function trigger') for t in selection.get('trigger_transcripts', [])]
            candidates += [(t, 'Broader review') for t in selection.get('transcripts', [])]
            for tx, role in candidates:
                if not tx:
                    continue
                item = dict(gene=tx.get('SYMBOL'), transcript=tx.get('Feature'), consequence=tx.get('Consequence'),
                            impact=tx.get('IMPACT'), hgvsc=tx.get('HGVSc'), hgvsp=tx.get('HGVSp'))
                same = next((v for v in transcripts if all(v.get(k) == value for k, value in item.items())), None)
                if same:
                    if role not in same['roles']:
                        same['roles'].append(role)
                else:
                    transcripts.append(dict(item, roles=[role]))
        tx = primary.get('gene_context', {}).get('transcript') or {}
        impact = primary.get('impact') or tx.get('IMPACT')
        high = [t for t in transcripts if t.get('impact') == 'HIGH']
        population = [e.get('population', {}) for e in row.get('evidence', [])]
        population += [s.get('population', {}) for s in selections]
        frequencies = [v for v in [*[s.get('population_af') for s in selections], *[p.get('combined_af') for p in population]] if isinstance(v, (int, float))]
        common = any(p.get('status') == 'excluded_common' for p in population)
        reasons = []
        for assertion in assertions:
            reason = 'ClinVar: ' + text(assertion.get('CLNSIG'))
            review = text(assertion.get('CLNREVSTAT'))
            reasons.append(reason + ('; ' + review if review else '; review status not available'))
        if conflict:
            reasons.append('Conflicting ClinVar assertions; pathogenicity is unresolved.')
        if uncertain:
            reasons.append('ClinVar includes uncertain significance; a predicted effect does not resolve that classification.')
        if benign:
            reasons.append('ClinVar includes benign / likely benign evidence; retained for supporting review.')
        if lof:
            effects = sorted({(t.get('gene') + ': ' if t.get('gene') else '') + text(t.get('consequence')) for t in transcripts if 'Loss-of-function trigger' in t['roles'] and t.get('consequence')})
            reasons.append('Predicted loss of function: ' + (', '.join(effects) or text(primary.get('consequence')) or 'effect not available') + '; HIGH impact on a transcript accepted by the gene-aware selection policy.')
        if panel_coding:
            effects = text(primary.get('consequence')) or 'protein-altering effect'
            reasons.append('Protein-altering variant in the reporting gene panel: ' + effects + '; VEP ' + (impact or 'impact not available') + ('.' if clean else '; pathogenicity is not established.'))
        if any(s.get('review_group') == 'consensus_review' for s in selections):
            effects = sorted({text(t.get('consequence')) for t in high if t.get('consequence')})
            criteria = ('VEP HIGH (' + ', '.join(effects) + ')' if high else 'ClinVar pathogenic / likely pathogenic assertion')
            reasons.append('Broader two-caller review: ' + criteria + '; ClairS-TO and DeepSomatic both passed the exact allele.')
        if frequencies:
            af = max(frequencies)
            reasons.append('Maximum reported population frequency: ' + f'{af*100:.6g}%' + ('.' if af <= .01 else '; above the 1% ordinary-review cutoff.'))
        else:
            reasons.append('Population frequency not available in the configured references.')
        if any(p.get('whitelist_rescue') for p in population):
            reasons.append('Exact allele is whitelisted; retained despite the ordinary population-frequency exclusion.' if frequencies and max(frequencies) > .01 else 'Exact allele is on the reporting whitelist.')
        if common:
            reasons.append('Common-population allele retained by broader consensus review, not the ordinary rare-variant criteria.')
        if any('panel_gene' in s for s in selections) and not any(s.get('panel_gene') for s in selections):
            reasons.append('Gene is outside the reporting panel.')
        if common or conflict or benign or (uncertain and plp):
            group, rank = 'supporting', 7 if common else 6
        elif clean:
            group, rank = 'findings', 0 if significance(sig, 'pathogenic') else 1
        elif lof:
            group, rank = 'findings', 2
        elif plp:
            group, rank = 'coding', 3
        elif high:
            group, rank = 'coding', 4
        elif panel_coding:
            group, rank = 'coding', 5
        else:
            group, rank = 'supporting', 6
        cosmic = [e.get('cosmic', {}) for e in row.get('evidence', [])]
        if any('targeted_evidence' in c for c in cosmic):
            matched = any(c.get('status') == 'matched' for c in cosmic)
            blood = any(c.get('blood_match') for c in cosmic)
            if matched:
                reasons.insert(0, 'Exact COSMIC allele match' + (' reported in blood/lymphoid tissue; primary review.' if blood else '; expanded review.'))
            for warning in sorted({w for c in cosmic for w in c.get('warnings', [])}):
                reasons.append('COSMIC join incomplete: ' + warning.replace('_', ' ') + '.')
            # Primary membership reflects tissue evidence, independently of ClinVar classification.
            if blood and not common:
                group = 'findings'
            elif not matched and any(c.get('status') == 'no_match' for c in cosmic) and group == 'findings':
                group = 'noncosmic_findings'
            else:
                group = 'coding'
            rank = 0 if blood else 1 if matched else 2 + rank
        if not native_reasons and not assertions:
            reasons.insert(0, 'Selected for review; original inclusion criteria are unavailable.')
        stars = max((4 if 'practice guideline' in text(a.get('CLNREVSTAT')) else 3 if 'expert panel' in text(a.get('CLNREVSTAT')) else 2 if 'multiple submitters' in text(a.get('CLNREVSTAT')) and 'no conflicts' in text(a.get('CLNREVSTAT')) else 1 if 'criteria provided' in text(a.get('CLNREVSTAT')) else 0 for a in assertions), default=0)
        row['presentation'] = dict(group=group, rank=rank, review_stars=stars, reasons=list(dict.fromkeys(reasons)),
                                   primary_selection=primary, selected_impact=impact, transcripts=transcripts,
                                   cosmic=cosmic_summary(row))
