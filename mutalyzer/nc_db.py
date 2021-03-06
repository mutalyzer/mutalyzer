from datetime import datetime

import mmap
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from mutalyzer.GenRecord import PList, Locus, Gene, Record
from mutalyzer.dbgb.models import Transcript, Reference
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

from mutalyzer.config import settings


def get_chromosome_ids(transcript_id):
    ids = []
    accession = transcript_id.split('.')[0]
    version = transcript_id.split('.')[1]

    references = Transcript.query. \
        filter_by(transcript_accession=accession). \
        filter_by(transcript_version=version).all()

    for transcript in references:
        reference = Reference.query. \
            filter_by(id=transcript.reference_id).first()
        ids.append('{}.{}'.format(reference.accession, reference.version))
    return ids


def get_entire_nc_record(record_id, geneName=None):

    # Get the accession
    accession = record_id.split('.')[0]
    version = record_id.split('.')[1]

    reference = _get_reference(accession, version)

    # If no reference present in the database we just return.
    if reference is None:
        return None

    if geneName is not None:
        db_transcripts = Transcript.query.\
            filter_by(reference_id=reference.id).\
            filter_by(gene=geneName).all()
    else:
        db_transcripts = _get_transcripts(reference, 1, reference.length)

    record = _get_mutalyzer_record(reference, db_transcripts)
    return record


def get_accession_version(reference):
    reference_list = str(reference).split('.')
    accession = version = None
    if len(reference_list) == 2:
        accession = reference_list[0]
        version = reference_list[1]
    if len(reference_list) == 1:
        accession = reference_list[0]

    return accession, version


def get_nc_record(record_id, parsed_description, output):
    """
    Get an NC record from the gbparser database and transform it into the
    Mutalyzer record format.

    :param record_id: HGVS format record description.
    :param parsed_description: HGVS description object.
    :param output: Mutalyzer special all purpose output object.
    :return: The record in mutalyzer format or None if not found.
    """
    if not _configuration_check(output):
        return None

    # Get the accession
    accession, version = get_accession_version(record_id)

    if version is None:
        versions = _get_versions(accession)
        if len(versions) > 0:
            output.addMessage(
                __file__, 2, 'ENOVERSION', 'Multiple versions found for %s. '
                                           'Please choose from: %s.'
                                           % (accession, ", ".join(versions)))
        return None

    reference = _get_reference(accession, version)

    # If no reference present in the database we just return.
    if reference is None:
        versions = _get_versions(accession)
        if len(versions) > 0:
            output.addMessage(
                __file__, 2, 'ENOVERSION', 'Version %s not found for %s '
                                           'accession in our database. You '
                                           'can choose instead from versions: '
                                           '%s.'
                                           % (version, accession,
                                              ", ".join(versions)))
        return None

    if parsed_description.RefType in ['g', '']:
        # Examples:
        # - NC_000001.11:g.111501128del
        # - NC_000001.11:62825del
        p_s, p_e = _get_description_boundary_positions(parsed_description)
    elif parsed_description.RefType in ['c', 'n'] or \
            parsed_description.Gene or parsed_description.AccNoTransVar:
        if parsed_description.Gene:
            # Example: 'NC_000001.11(OR4F5_v001):c.101del'
            transcripts = Transcript.query.\
                filter_by(reference_id=reference.id,
                          gene=parsed_description.Gene.GeneSymbol)\
                .all()
            if transcripts and len(transcripts) > 0:
                p_s, p_e = _boundaries(transcripts)
            else:
                return _record_with_genes_only(reference)
        elif parsed_description.AccNoTransVar:
            # Example: 'NC_000001.11(NM_001302680.1):c.62825del'
            accession = parsed_description.AccNoTransVar[0]
            version = parsed_description.AccNoTransVar[1]
            try:
                transcript = Transcript.query.\
                    filter_by(reference_id=reference.id,
                              transcript_accession=accession,
                              transcript_version=version).one()
            except NoResultFound:
                return _record_with_genes_only(reference)
            except MultipleResultsFound:
                output.addMessage(
                    __file__, 4, 'NCDBERROR', 'Possible NC DB error'
                                              ' - multiple entries.')
                return None
            else:
                p_s = transcript.transcript_start
                p_e = transcript.transcript_stop
        else:
            return _record_with_genes_only(reference)
    else:
        return _record_with_genes_only(reference)

    db_transcripts = _get_transcripts(reference, p_s, p_e)

    record = _get_mutalyzer_record(reference, db_transcripts)
    return record


def cds_position_list(mrna_position_list, cds_location):
    """
    Construct a list of coordinates that contains CDS start and stop and
    the internal splice sites.

    @arg mrna_position_list: mRNA positions/coordinates list
    @type mrna_position_list: list (integer)
    @arg cds_location: coding DNA positions/coordinates
    @type cds_location: list (integer)

    @return: CDS positions plus internal splice sites
    @rtype: list (integer)
    """
    i = 1
    ret = [cds_location[0]]

    while i <= len(mrna_position_list) - 1 \
            and cds_location[0] > mrna_position_list[i]:
        i += 2

    j = i
    while j <= len(mrna_position_list) \
            and cds_location[1] > mrna_position_list[j]:
        j += 2

    ret.extend(mrna_position_list[i:j])
    ret.append(cds_location[1])

    return ret


def _bare_record(reference):
    record = Record()
    # Populating the record with the generic information.
    record.source_id = reference.accession
    record.id = reference.accession
    record.source_accession = reference.accession
    record.source_version = reference.version
    record.organism = 'Homo sapiens'
    record.molType = 'g'
    return record


def _record_with_genes_only(reference):
    gene_names = Transcript.query.with_entities(Transcript.gene). \
        filter_by(reference_id=reference.id).all()
    record = _bare_record(reference)
    genes = []
    for gene_name in gene_names:
        gene = Gene(gene_name.gene)
        version = gene.newLocusTag()
        my_transcript = Locus(version)

        my_transcript.mRNA = PList()
        my_transcript.CDS = PList()

        gene.transcriptList.append(my_transcript)
        genes.append(gene)
    record.geneList = genes
    record.seq = Seq('', generic_dna)
    return record


def _get_mutalyzer_record(reference, db_transcripts):
    """
    Creates a Mutalyzer specific record from the transcript entries retrieved
    from the gbparser database.
    :param reference: A gbparser database reference entry.
    :param db_transcripts:A gbparser database list of transcript.
    :return: The Mutalyzer record.
    """
    record = _bare_record(reference)

    # Extracting the transcripts from the DB entries.
    transcripts = []
    for db_transcript in db_transcripts:
        transcript = {
            'gene': db_transcript.gene,
            'strand': db_transcript.strand,
            'transcript_start': db_transcript.transcript_start,
            'transcript_stop': db_transcript.transcript_stop,
            'transcript_product': db_transcript.transcript_product,
            'exons': [],
            'exons_start': db_transcript.exons_start,
            'exons_stop': db_transcript.exons_stop,
            'transcriptID': db_transcript.transcript_accession + '.' +
                            db_transcript.transcript_version,
        }
        if db_transcript.protein_accession is not None \
                and db_transcript.protein_version is not None:
            transcript['cds_start'] = db_transcript.cds_start
            transcript['cds_stop'] = db_transcript.cds_stop
            transcript['protein_product'] = db_transcript.protein_product
            transcript['proteinID'] = '%s.%s' %\
                                            (db_transcript.protein_accession,
                                             db_transcript.protein_version)
            transcript['linkMethod'] = 'ncbi'
        starts = map(int, db_transcript.exons_start.split(',')) \
            if db_transcript.exons_start else None
        stops = map(int, db_transcript.exons_stop.split(',')) \
            if db_transcript.exons_stop else None
        if (starts and stops) and (len(starts) == len(stops)):
            for start, stop in zip(starts, stops):
                exon = {'start': start,
                        'stop': stop}
                transcript['exons'].append(exon)
        transcripts.append(transcript)

    # Generating the actual record entries in the Mutalyzer format.
    gene_dict = {}
    for db_transcript in transcripts:
        if db_transcript['gene'] in gene_dict:
            gene = gene_dict[db_transcript['gene']]
        else:
            gene = Gene(db_transcript['gene'])

        if db_transcript['strand'] == '+':
            gene.orientation = 1
        if db_transcript['strand'] == '-':
            gene.orientation = -1

        transcript = Locus(gene.newLocusTag())

        transcript.mRNA = PList()
        transcript.mRNA.location = [db_transcript['transcript_start'],
                                       db_transcript['transcript_stop']]

        transcript.transcriptID = db_transcript['transcriptID']
        transcript.exon = PList()
        if db_transcript.get('exons') \
                and isinstance(db_transcript.get('exons'), list):
            exon_list = []
            for exon in db_transcript['exons']:
                exon_list.extend([exon['start'], exon['stop']])
            transcript.exon.positionList = exon_list
        else:
            transcript.exon.positionList = transcript.mRNA.location

        transcript.mRNA.positionList = transcript.exon.positionList
        transcript.mRNA.positionList.sort()

        if db_transcript.get('proteinID'):
            transcript.CDS = PList()
            transcript.CDS.location = [db_transcript['cds_start'],
                                          db_transcript['cds_stop']]

            transcript.CDS.positionList = cds_position_list(
                transcript.mRNA.positionList,
                transcript.CDS.location)

            transcript.proteinID = db_transcript['proteinID']

            transcript.transcriptProduct = db_transcript['transcript_product']
            transcript.proteinProduct = db_transcript['protein_product']
            transcript.linkMethod = 'ncbi'
            transcript.transcribe = True
            transcript.translate = True
        else:
            transcript.linkMethod = None
            transcript.transcribe = True
            transcript.translate = False
            transcript.locusTag = ''

        # transcript.molType = db_transcript['molType']

        gene.transcriptList.append(transcript)
        gene_dict[gene.name] = gene

    record.geneList = list(gene_dict.values())

    # Get the sequence.
    seq_path = settings.SEQ_PATH + reference.checksum_sequence + '.sequence'
    try:
        seq = Seq(_get_sequence_mmap(seq_path, 1, reference.length + 1),
                  generic_dna)
    except IOError:
        return None
    else:
        record.seq = seq

    return record


def _get_description_boundary_positions(description):
    """
    Determines the minimum and maximum positions that appear in the operations
    of an HGVS description.

    For 'NC_000001.11:g.[100T>C;750A>G;2000C>T]' it should return 100 and 2000.

    :param description: Parsed HGVS description.
    :return: Minimum and maximum positions that appear in the description.
    """
    if description.SingleAlleleVarSet:
        variants = [v.RawVar for v in description.SingleAlleleVarSet]
    else:
        variants = [description.RawVar]

    positions = set()
    for variant in variants:
        first_location = last_location = variant.StartLoc.PtLoc
        if variant.EndLoc:
            last_location = variant.EndLoc.PtLoc
        positions.add(int(first_location.Main))
        positions.add(int(last_location.Main))

    return min(positions), max(positions)


def _configuration_check(output):
    # Check if the required variables were set into settings.
    try:
        settings.SEQ_PATH
        settings.DATABASE_GB_URI
    except (NameError, AttributeError):
        output.addMessage(
            __file__, 2, 'NCSETTINGS', 'Chromosomal DB settings not set.')
        return False

    # Check if the sequence folder exists.
    if not os.path.isdir(settings.SEQ_PATH):
        output.addMessage(
            __file__, 2, 'NCSEQDIR', 'Sequence directory does not exist.')
        return False

    return True


def _get_db_boundaries_positions(reference, position_start, position_end):
    """
    Providing two positions on a chromosome reference this method will query
    the gbparser database and will return the extremity positions of the
    transcript entries found that contain the input positions.

    Example 1:
                      p_s                        p_e
                       |                          |
            |------------------------------------------------|
                                                |---------------------|
           t_1                                                          t_2

    Example 2:
                      p_s                        p_e
                       |                          |
            |------------------|
                                             |---------------------|
           t_1                                                    t_2

    Example 2:
                      p_s                        p_e
                       |                          |
                           |------------------|
                                          |---------------------|
                      t_1                                      t_2

    return: t_1, t_2

    :param reference: Database reference entry.
    :param position_start:
    :param position_end:
    :return:
    """
    p_s = position_start
    p_e = position_end

    transcripts = Transcript.query.\
        filter_by(reference_id=reference.id). \
        filter(((Transcript.transcript_start <= p_s) &
                (Transcript.transcript_stop >= p_e)) |
               ((Transcript.transcript_start >= p_s) &
                (Transcript.transcript_stop <= p_e)) |
               ((Transcript.transcript_start <= p_s) &
                (Transcript.transcript_stop >= p_s)) |
               ((Transcript.transcript_start <= p_e) &
                (Transcript.transcript_stop >= p_e)))\
        .all()

    if transcripts and len(transcripts) > 0:
        return _boundaries(transcripts)
    else:
        return p_s, p_e


def _boundaries(transcripts):
    positions = set()

    for transcript in transcripts:
        positions.add(transcript.transcript_start)
        positions.add(transcript.transcript_stop)

    p_s = min(positions)
    p_e = max(positions)

    return p_s, p_e


def _get_reference(accession, version):
    """
    Retrieves the database reference entry for the user provided accession and
    version.
    :param accession: The accession from the user provided description.
    :param version: The version number from the user provided description.
    :return: Database reference entry.
    """
    reference = Reference.query.\
        filter_by(accession=str(accession),
                  version=str(version))\
        .order_by(Reference.id.asc()) \
        .first()
    return reference


def _get_versions(accession):
    """
    Retrieves a list with all the versions corresponding to an accession that
    is present in the database.
    :param accession: The accession for which to look for the versions.
    :return: List with the versions for the provided accession.
    """
    references = Reference.query.\
        filter_by(accession=str(accession))\
        .order_by(Reference.id.asc()) \
        .all()
    versions = []
    for reference in references:
        versions.append(reference.version)
    return versions


def _get_transcripts(reference, position_start, position_end):
    """
    Retrieves the transcripts information from the database for the provided
    reference that are between the provided start and end positions to which
    5000 is subtracted and added, respectively.
    :param reference:
    :param position_start:
    :param position_end:
    :return:
    """
    if position_start > 5000:
        p_s = position_start - 5000
    else:
        p_s = 1

    if position_end < reference.length - 5000:
        p_e = position_end + 5000
    else:
        p_e = reference.length

    p_s, p_e = _get_db_boundaries_positions(reference, p_s, p_e)

    transcripts = Transcript.query.filter_by(reference_id=reference.id). \
        filter(((Transcript.transcript_start <= p_s) &
                (Transcript.transcript_stop >= p_e)) |
               ((Transcript.transcript_start >= p_s) &
                (Transcript.transcript_stop <= p_e)) |
               ((Transcript.transcript_start <= p_s) &
                (Transcript.transcript_stop >= p_s)) |
               ((Transcript.transcript_start <= p_e) &
                (Transcript.transcript_stop >= p_e))).all()

    return transcripts


def _get_sequence_mmap(file_path, start, end):
    """
    Sequence retrieval.
    :param file_path: Path towards the sequence file.
    :param start: Start position.
    :param end: End position
    :return: The sequence.
    """
    with open(file_path, "r+b") as f:
        # memory-map the file, size 0 means whole file
        mm = mmap.mmap(f.fileno(), 0)
        return mm[start - 1:end]
