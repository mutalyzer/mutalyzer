import json
from .util import get_start, get_end


def get_mol_type(reference):
    if reference['model'].get('qualifiers'):
        if reference['model']['qualifiers'].get('mol_type'):
            return reference['model']['qualifiers']['mol_type']


def get_transcripts_ids(reference_model):
    transcript_ids = set()
    for feature in reference_model:
        if feature['type'] == 'gene':
            for sub_feature in feature['sub_features']:
                if sub_feature['type'] == 'mRNA':
                    transcript_ids.add(sub_feature['id'].split('-')[1])
    return list(transcript_ids)


def get_selector_model(reference_model, mol_type, selector_id=None):
    if 'DNA' in mol_type.upper():
            exons, cds = get_exon_cds_genomic(
                selector_id, reference_model['model'])
    elif mol_type == 'mRNA':
        exons, cds = get_exon_cds_for_mrna_reference(
            reference_model['model'])
    cds = sorted(cds)
    if len(cds) >= 2:
        cds = sorted([cds[0], cds[-1]])
    return {'exons': sorted(exons), 'cds': cds}


def get_exons_cds(feature):
    """
    Get the exons and the cds for a provided feature.
    """
    exons = []
    cds = []
    if feature.get('features'):
        for sub_feature in feature['features']:
            if sub_feature['type'] == 'exon':
                exons.append((get_start(sub_feature), get_end(sub_feature)))
            elif sub_feature['type'] == 'CDS':
                cds.extend([get_start(sub_feature), get_end(sub_feature)])
    cds = sorted(cds)
    if len(cds) >= 2:
        cds = sorted([cds[0], cds[-1]])
    return sorted(exons), cds


def get_exon_cds_genomic_ncbi(selector_id, reference_model):
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
                            return get_exons_cds(sub_feature)
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
                            return get_exons_cds(sub_feature)
    return [], []


def get_exon_cds_genomic(selector_id, reference_model):
    for feature in reference_model['features']:
        if feature['type'] == 'gene' and feature.get('features'):
            for sub_feature in feature['features']:
                if sub_feature['id'] == selector_id:
                    return get_exons_cds(sub_feature)
    return [], []


def get_all_selectors_exon_cds(reference_model):
    output = []
    if reference_model.get('features'):
        for feature in reference_model['features']:
            if feature['type'] == 'gene' and feature.get('features'):
                rna_index = 1
                for sub_feature in feature['features']:
                    gene_id = '{}_v{:03}'.format(feature['id'], rna_index)
                    rna_id = sub_feature['id']
                    exons, cds = get_exons_cds(sub_feature)
                    if exons:
                        output.append({'exons': exons,
                                       'cds': cds,
                                       'id1': gene_id,
                                       'id2': rna_id})
                    rna_index += 1
    return output


def get_exon_cds_for_mrna_reference(reference_model):
    if reference_model.get('features') and reference_model['features'][0] and \
            reference_model['features'][0].get('features'):
        return get_exons_cds(reference_model['features'][0].get('features')[0])
    return [], []
