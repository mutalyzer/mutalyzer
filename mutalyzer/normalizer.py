# from extractor import describe_dna
# from mutalyzer_mutator import mutate
# from mutalyzer_mutator.util import reverse_complement
#
# import mutalyzer.infos as infos
#
# from .converter.extras import convert_reference_model
# from .converter.to_hgvs_coordinates import to_hgvs_locations
# from .converter.variants_de_to_hgvs import de_to_hgvs
from .description import Description

# from .description_model import get_reference_id
# from .reference import (
#     get_chromosome_accession,
#     get_coordinate_system_from_selector_id,
#     get_internal_selector_model,
#     get_reference_mol_type,
#     get_selectors_ids,
#     retrieve_reference,
# )
#
#
# def get_chromosomal_description(d):
#     # TODO: Add tests.
#
#     if d.errors or not d.references or d.only_variants:
#         return
#     ref_id = get_reference_id(d.corrected_model)
#     if (
#         ref_id
#         and get_reference_mol_type(d.references[ref_id]) == "mRNA"
#         and d.corrected_model["coordinate_system"] == "c"
#         and (ref_id.startswith("NM_") or ref_id.startswith("XM_"))
#     ):
#         d.add_info(infos.mrna_genomic_tip())
#     elif (
#         ref_id
#         and get_reference_mol_type(d.references[ref_id]) == "mRNA"
#         and d.corrected_model["coordinate_system"] == "r"
#         and (ref_id.startswith("NM_") or ref_id.startswith("XM_"))
#     ):
#         return
#
#     chromosome_accessions = get_chromosome_accession(ref_id, d.references["reference"])
#     if not chromosome_accessions:
#         return
#
#     chromosomal_descriptions = []
#     for assembly, chromosome_accession in chromosome_accessions:
#         chromosome_model = retrieve_reference(chromosome_accession)[0]
#         if chromosome_model:
#             selector_ids = get_selectors_ids(chromosome_model["annotations"], "c")
#             ref_id_accession = ref_id.split(".")[0]
#             matched_selector_ids = [i for i in selector_ids if i.startswith(ref_id_accession)]
#
#             if not matched_selector_ids:
#                 d.add_info(infos.no_selector(chromosome_accession, ref_id))
#                 continue
#             if ref_id not in matched_selector_ids:
#                 d.add_info(infos.other_versions(chromosome_accession, ref_id, matched_selector_ids))
#                 continue
#
#             selector_id = ref_id
#
#             from_reference_model = convert_reference_model(d.references["reference"], selector_id, "transcript")
#
#             ref_seq_from = from_reference_model["sequence"]["seq"]
#             seqs = d.get_sequences()
#             seqs["reference"] = from_reference_model["sequence"]["seq"]
#             seqs[selector_id] = from_reference_model["sequence"]["seq"]
#             obs_seq = mutate(seqs, d.delins_model["variants"])
#
#             to_reference_model = convert_reference_model(
#                 chromosome_model, selector_id, "transcript"
#             )
#
#             ref_seq_to = to_reference_model["sequence"]["seq"]
#
#             to_inverted = get_internal_selector_model(
#                 chromosome_model["annotations"], selector_id, True
#             ).get("inverted")
#             if to_inverted:
#                 obs_seq = reverse_complement(obs_seq)
#             variants = de_to_hgvs(
#                 describe_dna(ref_seq_to, obs_seq),
#                 {"reference": ref_seq_to, "observed": obs_seq},
#             )
#
#             if (not to_inverted and ref_seq_to != ref_seq_from) or (
#                 to_inverted and reverse_complement(ref_seq_to) != ref_seq_from
#             ):
#                 d.add_info(
#                     infos.mrna_genomic_difference(ref_id, chromosome_accession)
#                 )
#             else:
#                 variants_model = to_hgvs_locations(
#                     {
#                         "reference": {
#                             "id": chromosome_accession,
#                             "selector": {"id": selector_id},
#                         },
#                         "coordinate_system": "i",
#                         "variants": variants,
#                     },
#                     {
#                         "reference": to_reference_model,
#                         chromosome_model["annotations"]["id"]: to_reference_model,
#                     },
#                     get_coordinate_system_from_selector_id(
#                         chromosome_model, selector_id
#                     ),
#                     selector_id,
#                     True,
#                 )
#                 chr_d = Description(description_model=variants_model)
#                 chr_d.to_delins()
#                 chr_d.mutate()
#                 chr_d.extract()
#                 chr_d.construct_de_hgvs_internal_indexing_model()
#                 chr_d.construct_de_hgvs_coordinates_model()
#                 chr_d.construct_normalized_description()
#
#                 chromosomal_descriptions.append(
#                     {
#                         "assembly": assembly,
#                         "description": chr_d.normalized_description,
#                     }
#                 )
#
#     if chromosomal_descriptions:
#         d.chromosomal_descriptions = chromosomal_descriptions
#


def normalize(description, only_variants=False, sequence=None):
    d = Description(
        description=description,
        only_variants=only_variants,
        sequence=sequence,
    )
    d.normalize()

    # get_chromosomal_description(d)

    output = d.output()

    return output
