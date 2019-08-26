import json


def get_transcripts_ids(reference_model):
    transcript_ids = set()
    for feature in reference_model:
        if feature['type'] == 'gene':
            for sub_feature in feature['sub_features']:
                if sub_feature['type'] == 'mRNA':
                    transcript_ids.add(sub_feature['id'].split('-')[1])
    return list(transcript_ids)
