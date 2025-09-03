from os.path import commonprefix

from algebra import Variant
from mutalyzer_mutator.util import reverse_complement

from mutalyzer.util import (
    create_exact_point_model,
    create_exact_range_model,
    get_end,
    get_inserted_sequence,
    get_start,
)

from .util import trim


def get_dominators(graph):
    successors = {}
    all_nodes = set()
    tail_nodes = set()
    head_nodes = set()

    for edge in graph.edges():
        head, tail = edge["head"], edge["tail"]
        all_nodes.add(head)
        all_nodes.add(tail)
        head_nodes.add(head)
        tail_nodes.add(tail)

        if head not in successors:
            successors[head] = []
        successors[head].append(tail)

    sink = next(iter(tail_nodes - head_nodes), None)
    source = next(iter(head_nodes - tail_nodes), None)

    if sink is None or source is None:
        return set()

    source = (head_nodes - tail_nodes).pop()
    sink = (tail_nodes - head_nodes).pop()

    dominators = {source: {source}}

    for node in all_nodes:
        if node != source:
            dominators[node] = all_nodes.copy()

    changed = True
    while changed:
        changed = False
        for node in all_nodes:
            if node == source:
                continue

            predecessors = [n for n in all_nodes if n in successors and node in successors[n]]

            if predecessors:
                new_dominators = set.intersection(*[dominators[pred] for pred in predecessors])
                new_dominators.add(node)

                if new_dominators != dominators[node]:
                    dominators[node] = new_dominators
                    changed = True

    return dominators[sink] - {source, sink}


def graph_to_dot(graph, reference, labels=True, dominators=True, complexity_limit=1000):
    width = ".8" if labels else "1"

    dot_lines = [
        "digraph {",
        "rankdir=LR",
        "edge[fontname=monospace]",
        f'node[fixedsize=true,fontname=serif,shape=circle,width={width}]'
    ]

    nodes = {}
    head_nodes = set()
    tail_nodes = set()
    node_index = 0
    edge_index = 0

    for edge in graph.edges():
        edge_index += 1

        if edge_index > 200 or edge_index * len(nodes) > complexity_limit:
            return f"// Graph too complex. Stopped at {edge_index} edges and {len(nodes)} nodes."

        head, tail, variant, count = edge["head"], edge["tail"], edge["variant"], edge["count"]

        tail_nodes.add(tail)
        head_nodes.add(head)

        if tail not in nodes:
            nodes[tail] = f"s{node_index}"
            node_index += 1
        if head not in nodes:
            nodes[head] = f"s{node_index}"
            node_index += 1

        head_id, tail_id = nodes[head], nodes[tail]

        if variant:
            label = to_hgvs(variant, reference)
            if count > 1:
                dot_lines.append(f'  {head_id} -> {tail_id} [label="{label} x {count}",penwidth=2]')
            else:
                dot_lines.append(f'  {head_id} -> {tail_id} [label="{label}"]')
        else:
            dot_lines.append(f'  {head_id} -> {tail_id} [label="&lambda;",style=dashed]')

    sink_node = next(iter(tail_nodes - head_nodes), None)
    source_node = next(iter(head_nodes - tail_nodes), None)

    if sink_node:
        dot_lines.append(f'{nodes[sink_node]}[fillcolor=aliceblue,style=filled,peripheries=2,penwidth=2]')

    if source_node:
        dot_lines.extend([
            f'{nodes[source_node]}[fillcolor=aliceblue,style=filled,penwidth=2]',
            'i[label="",shape=point,width=.1]',
            f'i->{nodes[source_node]}',
        ])

    if dominators:
        for node in get_dominators(graph):
            dot_lines.append(f'{nodes[node]}[fillcolor=aliceblue,style=filled,penwidth=2]')
    dot_lines.append("}")
    return "\n".join(dot_lines)


def trim_for_reverse(lhs, rhs):
    """Find the lengths of the common prefix and common suffix between
    two sequences."""
    idx = len(commonprefix([lhs[::-1], rhs[::-1]]))
    return len(commonprefix([lhs[:len(lhs)-idx], rhs[:len(rhs)-idx]])), idx


def to_hgvs_dict(variants, ref_seq, forward_strand=True):
    """Algebra based experimental version of HGVS serialization with support for
    tandem repeats and complex variants."""
    def var_dict(var_type, start, end=None, inserted=None, repeat_number=None, del_seq=None):
        output = {
            "location": to_hgvs_position(start, end),
            "type": var_type,
            "source": "reference",
        }
        if isinstance(inserted, list):
            output["inserted"] = inserted
        else:
            if inserted:
                output["inserted"] = [{"sequence": inserted, "source": "description"}]
            if inserted and repeat_number:
                output["inserted"][0]["repeat_number"] = {"type": "point", "value": repeat_number}
        if del_seq:
            output["deleted"] = [{"sequence": del_seq, "source": "description"}]
        return output

    def repeats(word):
        length = 0
        idx = 1
        lps = [0] * len(word)
        while idx < len(word):
            if word[idx] == word[length]:
                length += 1
                lps[idx] = length
                idx += 1
            elif length != 0:
                length = lps[length - 1]
            else:
                lps[idx] = 0
                idx += 1

        pattern = len(word) - length
        if pattern == 0:
            return "", 0, 0
        return word[:pattern], len(word) // pattern, len(word) % pattern

    def to_hgvs_position(start, end=None):
        if end is None or end - start == 1:
            return create_exact_point_model(start + 1)
        if start == end:
            return create_exact_range_model(start, start + 1)
        return create_exact_range_model(start + 1, end)

    def other(variant):
        if variant.end - variant.start == 0:
            if not variant.sequence:
                return "="
            # print("insertion sequence", variant.sequence)
            return var_dict("insertion", variant.start, variant.start, variant.sequence)
            # return f"{variant.start}_{variant.start + 1}ins{variant.sequence}"

        deleted = ""
        substitution = ref_seq[variant.start:variant.end]

        if variant.end - variant.start == 1:
            if not variant.sequence:
                return var_dict("deletion", variant.start)
                # return f"{variant.start + 1}del{deleted}"
            if len(variant.sequence) == 1:
                return var_dict("substitution", variant.start, del_seq=substitution, inserted=variant.sequence)
                # return f"{variant.start + 1}{substitution}>{variant.sequence}"
            return var_dict("deletion_insertion", variant.start, del_seq=deleted, inserted=variant.sequence)
            # return f"{variant.start + 1}del{deleted}ins{variant.sequence}"

        if not variant.sequence:
            return var_dict("deletion", variant.start, variant.end, del_seq=deleted)
            # return f"{variant.start + 1}_{variant.end}del{deleted}"

        return var_dict("deletion_insertion", variant.start, variant.end, del_seq=deleted, inserted=variant.sequence)
        # return f"{variant.start + 1}_{variant.end}del{deleted}ins{variant.sequence}"

    def hgvs(variant):
        inserted_unit, inserted_number, inserted_remainder = repeats(variant.sequence)
        deleted = ref_seq[variant.start:variant.end]
        deleted_unit, deleted_number, deleted_remainder = repeats(deleted)

        # Select a non-minimal repeat unit if reference and observed are
        # in agreement.
        diff = len(inserted_unit) - len(deleted_unit)
        if diff < 0 and deleted_unit == variant.sequence[:len(inserted_unit) - diff]:
            inserted_unit = deleted_unit
            inserted_number = 1
            inserted_remainder = deleted_remainder
        elif diff > 0 and inserted_unit == deleted[:len(deleted_unit) + diff]:
            deleted_unit = inserted_unit
            deleted_number = 1
            deleted_remainder = inserted_remainder

        # print(f"{to_hgvs_position(variant.start, variant.end - alt_deleted_remainder)}{inserted_unit}[{inserted_number}]")

        # Repeat structure
        if deleted_unit == inserted_unit:
            if deleted_number == inserted_number:
                raise ValueError("empty variant")

            # Duplication
            if deleted_number == 1 and inserted_number == 2:
                if forward_strand:
                    return var_dict(
                        "duplication",
                        variant.start + inserted_remainder,
                        variant.start + inserted_remainder + len(inserted_unit),
                    )
                    # return f"{to_hgvs_position(variant.start + inserted_remainder, variant.start + inserted_remainder + len(inserted_unit))}dup"
                return var_dict(
                    "duplication",
                    variant.start,
                    variant.start + len(inserted_unit),
                )
            # shift 3'
            assert deleted_remainder == inserted_remainder
            if forward_strand:
                inserted_unit = variant.sequence[inserted_remainder:inserted_remainder + len(inserted_unit)]
                r_d = var_dict("repeat", variant.start + deleted_remainder, variant.end, inserted_unit, inserted_number)
                if r_d.get("location", {}).get("start"):
                    r_d["location"]["start"]["shift"] = deleted_remainder
                if r_d.get("location", {}).get("end"):
                    r_d["location"]["end"]["shift"] = deleted_remainder
                return r_d
            inserted_unit = variant.sequence[:len(inserted_unit)]
            return var_dict("repeat", variant.start, variant.end - deleted_remainder, inserted_unit, inserted_number)
            # return f"{to_hgvs_position(variant.start + deleted_remainder, variant.end)}{inserted_unit}[{inserted_number}]"

        # Prefix and suffix trimming
        if forward_strand:
            start, end = trim(deleted, variant.sequence)
        else:
            start, end = trim_for_reverse(deleted, variant.sequence)
        trimmed = Variant(variant.start + start, variant.end - end, variant.sequence[start:len(variant.sequence) - end])

        # Inversion
        if len(trimmed.sequence) > 1 and trimmed.sequence == reverse_complement(ref_seq[trimmed.start:trimmed.end]):
            return var_dict("inversion", trimmed.start, trimmed.end)
            # return f"{to_hgvs_position(trimmed.start, trimmed.end)}inv"

        # Deletion/insertion with repeated insertion
        inserted_unit, inserted_number, inserted_remainder = repeats(trimmed.sequence)
        if inserted_number > 1:
            suffix = [{"sequence": inserted_unit, "source": "description", "repeat_number": {"type": "point", "value": inserted_number}}]
            # suffix = f"{inserted_unit}[{inserted_number}]"
            if inserted_remainder:
                suffix = [suffix[0], {"sequence": inserted_unit[:inserted_remainder], "source": "description"}]
                # suffix = f"[{suffix};{inserted_unit[:inserted_remainder]}]"

            if trimmed.start == trimmed.end:
                if forward_strand:
                    return var_dict("insertion", trimmed.start, trimmed.start, suffix)
                # return f"{to_hgvs_position(trimmed.start, trimmed.end)}ins{suffix}"
            return var_dict("deletion_insertion", trimmed.start, trimmed.end, suffix)
            # return f"{to_hgvs_position(trimmed.start, trimmed.end)}delins{suffix}"

        # All other variants
        return other(trimmed)
        # return trimmed.to_hgvs(ref_seq)

    if not variants:
        return []

    if len(variants) == 1:
        return [hgvs(variants[0])]

    return [hgvs(variant) for variant in variants]


def algebra_variants(variants_delins, sequences):
    variants_algebra = []
    for variant in variants_delins:
        variants_algebra.append(
            Variant(get_start(variant), get_end(variant), get_inserted_sequence(variant, sequences))
        )
    return variants_algebra



def to_hgvs(variant, reference=None, only_substitutions=True):
    """
    Taken from the algebra.
    """

    if variant.end - variant.start == 0:
        if not variant.sequence:
            return "="
        return f"{variant.start}_{variant.start + 1}ins{variant.sequence}"

    deleted = ""
    substitution = ""
    if reference is not None:
        if not only_substitutions:
            deleted = reference[variant.start:variant.end]
        substitution = reference[variant.start:variant.end]

    if variant.end - variant.start == 1:
        if not variant.sequence:
            return f"{variant.start + 1}del{deleted}"
        if len(variant.sequence) == 1:
            return f"{variant.start + 1}{substitution}>{variant.sequence}"
        return f"{variant.start + 1}del{deleted}ins{variant.sequence}"

    if not variant.sequence:
        return f"{variant.start + 1}_{variant.end}del{deleted}"

    return f"{variant.start + 1}_{variant.end}del{deleted}ins{variant.sequence}"


def to_spdi(variant, reference_id=""):
    """
    Taken from the algebra.
    """
    return f"{reference_id}:{variant.start}:{variant.end - variant.start}:"f"{variant.sequence}"


def algebra_variant_to_delins(variant):
    delins_variant = {
        "type": "deletion_insertion",
        "source": "reference",
        "location": {
            "type": "range",
            "start": {"type": "point", "position": variant.start},
            "end": {"type": "point", "position": variant.end},
        },
        "deleted": [],
        "inserted": [],
    }
    if variant.sequence:
        delins_variant["inserted"].append(
            {"sequence": variant.sequence, "source": "description"}
        )

    return delins_variant


def delins_to_algebra_variant(v, sequences):
    return Variant(get_start(v), get_end(v), get_inserted_sequence(v, sequences))


def delins_to_algebra(variants, sequences):
    return [delins_to_algebra_variant(v, sequences) for v in variants]
