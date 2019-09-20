import json
from .util import get_start, get_end


def get_mol_type(reference):
    if reference['source'] == 'lrg':
        return 'genomic DNA'
    for feature in reference['model']['features']:
        if feature['type'] == 'region':
            return feature['qualifiers'].get('mol_type')


def get_transcripts_ids(reference_model):
    transcript_ids = set()
    for feature in reference_model:
        if feature['type'] == 'gene':
            for sub_feature in feature['sub_features']:
                if sub_feature['type'] == 'mRNA':
                    transcript_ids.add(sub_feature['id'].split('-')[1])
    return list(transcript_ids)


def get_selector_model(reference_model, mol_type, selector_id=None):
    if mol_type == 'genomic DNA':
        if reference_model['source'] == 'ncbi':
            exons, cds = get_exon_cds_genomic_ncbi(
                selector_id, reference_model['model'])
        elif reference_model['source'] == 'lrg':
            exons, cds = get_exon_cds_genomic_lrg(
                selector_id, reference_model['model'])
    elif mol_type == 'mRNA':
        exons, cds = get_exon_cds_for_mrna_reference(
            reference_model['model'])
    cds = sorted(cds)
    if len(cds) >= 2:
        cds = sorted([cds[0], cds[-1]])
    return {'exons': sorted(exons), 'cds': cds}


def get_exon_cds_genomic_ncbi(selector_id, reference_model):
    exons = []
    cds = []
    if '_v' in selector_id:
        gene_id = selector_id.split('_v')[0]
        transcript_number = int(selector_id.split('_v')[1])
        for feature in reference_model['features']:
            if feature['type'] == 'gene' and feature.get('features') and \
                    '-' in feature['id'] and \
                    feature['id'].split('gene-')[1] == gene_id:
                rna_id = 1
                for sub_feature in feature['features']:
                    if 'RNA' in sub_feature['id'].upper():
                        if rna_id == transcript_number:
                            for part in sub_feature['features']:
                                if part['type'] == 'exon':
                                    exons.append(
                                        (get_start(part), get_end(part)))
                                elif part['type'] == 'CDS':
                                    cds.extend(
                                        [get_start(part), get_end(part)])
                        rna_id += 1
    else:
        for feature in reference_model['features']:
            if feature['type'] == 'gene':
                if feature.get('features'):
                    for sub_feature in feature['features']:
                        if 'RNA' in sub_feature['type'].upper() and \
                                '-' in sub_feature['id'] and \
                                selector_id == \
                                sub_feature['id'].split('-')[1]:
                            for part in sub_feature['features']:
                                if part['type'] == 'exon':
                                    exons.append(
                                        (get_start(part), get_end(part)))
                                elif part['type'] == 'CDS':
                                    cds.extend(
                                        [get_start(part), get_end(part)])
    return exons, cds


def get_exon_cds_genomic_lrg(selector_id, reference_model):
    exons = []
    cds = []
    for feature in reference_model['features']:
        if feature['type'] == 'gene' and feature.get('features'):
            for sub_feature in feature['features']:
                if sub_feature['id'] == selector_id:
                    for part in sub_feature['features']:
                        if part['type'] == 'exon':
                            exons.append((get_start(part), get_end(part)))
                        elif part['type'] == 'cds':
                            cds.extend([get_start(part), get_end(part)])
    return exons, cds


def get_all_exon_cds_for_genomic(reference_model):
    output = []
    for feature in reference_model['features']:
        if feature['type'] == 'gene' and feature.get('features') and \
                '-' in feature['id']:
            rna_index = 1
            for sub_feature in feature['features']:
                gene_id = '{}_v{:03}'.format(
                    feature['id'].split('gene-')[1], rna_index)
                if 'rna' in sub_feature['id']:
                    rna_id = sub_feature['id'].split('-')[1]
                    exons = []
                    cds = []
                    for part in sub_feature['features']:
                        if part['type'] == 'exon':
                            exons.append((get_start(part), get_end(part)))
                        elif part['type'] == 'CDS':
                            cds.extend([get_start(part), get_end(part)])
                    if len(cds) >= 2:
                        cds = sorted([cds[0], cds[-1]])
                    else:
                        cds = []
                    output.append({'exons': exons,
                                   'cds': cds,
                                   'id1': gene_id,
                                   'id2': rna_id})
                    rna_index += 1
    return output


def get_exon_cds_for_mrna_reference(reference_model):
    exons = []
    cds = []
    for feature in reference_model['features']:
        if feature['type'] == 'gene' and feature.get('features'):
            for sub_feature in feature['features']:
                if sub_feature['type'] == 'CDS':
                    cds.append(sub_feature['location']['start']['position'])
                    cds.append(sub_feature['location']['end']['position'])
                elif sub_feature['type'] == 'exon':
                    exons.append((sub_feature['location']['start']['position'],
                                  sub_feature['location']['end']['position']))
    return exons, cds
