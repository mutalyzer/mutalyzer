from os.path import commonprefix

from algebra import Variant
from mutalyzer_mutator.util import reverse_complement

from mutalyzer.converter.extras import g_to_cn
from mutalyzer.util import (
    create_exact_point_model,
    create_exact_range_model,
    get_end,
    get_inserted_sequence,
    get_start,
)

from .util import trim


def algebra_variant_to_name_model(variant):
    def _position_to_hgvs():
        if variant.end - variant.start == 1:
            return {"type": "point", "position": variant.start + 1}
        if variant.start == variant.end:
            return {
                "type": "range",
                "start": {"type": "point", "position": variant.start},
                "end": {"type": "point", "position": variant.start + 1},
            }
        return {
            "type": "range",
            "start": {"type": "point", "position": variant.start + 1},
            "end": {"type": "point", "position": variant.end},
        }

    delins_variant = {
        "source": "reference",
        "location": _position_to_hgvs(),
        "deleted": [],
        "inserted": [],
    }
    if variant.sequence:
        if variant.start == variant.end:
            delins_variant["type"] = "insertion"
        else:
            delins_variant["type"] = "deletion_insertion"
        delins_variant["inserted"].append(
            {"sequence": variant.sequence, "source": "description"}
        )
    else:
        delins_variant["type"] = "deletion"
    return delins_variant


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


def variant_label(variant, count, reference, selector):
    """
    Generate HGVS label(s) for variant(s).

    Args:
        variant: Variant object with start, end, and sequence attributes.
        count: Number of variants to display.
        reference: Reference sequence.
        selector: Controls output order (reversed if last element is True).

    Returns:
        String with HGVS notation for the variant(s).
    """
    label = to_hgvs(variant, reference, selector)

    if count <= 1:
        return label

    last_sequence = _rotate_sequence(variant.sequence, count - 1)
    variant_last = Variant(variant.start + (count - 1), variant.end + (count - 1), last_sequence)
    label_end = to_hgvs(variant_last, reference, selector)

    if count == 3:
        middle_sequence = _rotate_sequence(variant.sequence, count - 2)
        variant_middle = Variant(variant.start + (count - 2), variant.end + (count - 2), middle_sequence)
        middle_content = to_hgvs(variant_middle, reference, selector)
    elif count > 3:
        middle_content = f"(+{count - 2} more)"
    else:  # count == 2
        middle_content = ""

    parts = [label, label_end]
    if middle_content:
        parts = [label, middle_content, label_end]

    if selector and len(selector) > 0 and selector[-1]:
        parts = parts[::-1]

    return "\n".join(parts)


def _rotate_sequence(sequence, offset):
    """Rotate sequence by offset positions."""
    if not sequence:
        return ""
    offset = offset % len(sequence)
    return sequence[offset:] + sequence[:offset]


def _assign_node_ids(all_nodes, source_node, sink_node):
    """
    Assign DOT node IDs with source as s0 and sink as highest number.
    """
    nodes = {}
    node_info = {}
    node_index = 0
    sorted_nodes = sorted(all_nodes)

    for node in sorted_nodes:
        if node == source_node:
            nodes[node] = "s0"
        elif node == sink_node:
            # Reserve the last ID for sink
            continue
        else:
            # Skip s0 since it's reserved for source
            if source_node is not None:
                nodes[node] = f"s{node_index + 1}"
                node_index += 1
            else:
                nodes[node] = f"s{node_index}"
                node_index += 1

        node_info[node] = f"{node}"

    # Assign sink node the highest number
    if sink_node is not None:
        if source_node is not None:
            nodes[sink_node] = f"s{node_index + 1}"
        else:
            nodes[sink_node] = f"s{node_index}"
        node_info[sink_node] = f"{sink_node}"

    return nodes, node_info


def _collect_nodes(graph, selector, edges_limit):
    """
    Collect all nodes and identify special nodes (source/sink).
    """
    all_nodes = set()
    head_nodes = set()
    tail_nodes = set()
    edge_count = 0
    limit_exceeded = False

    for edge in graph.edges():
        edge_count += 1

        head, tail = edge["head"], edge["tail"]

        if selector and selector[-1] is True:
            head, tail = tail, head

        all_nodes.add(head)
        all_nodes.add(tail)
        head_nodes.add(head)
        tail_nodes.add(tail)

        if edge_count >= edges_limit:
            limit_exceeded = True
            break

    return all_nodes, head_nodes, tail_nodes, edge_count, limit_exceeded


def _create_edge_dot(head, tail, variant, count, reference, selector, nodes):
    """
    Create DOT notation for a single edge.
    """
    if variant:
        label = variant_label(variant, count, reference, selector)
        pen_width = "1" if count == 1 else "2"
        tooltip = f"Variant: {variant}\\nCount: {count}"

        return (
            f'  {nodes[head]} -> {nodes[tail]} '
            f'[label="{label}",penwidth={pen_width},'
            f'tooltip="{tooltip}",edgetooltip="{tooltip}",labeltooltip="{tooltip}"]'
        )
    else:
        empty_tooltip = "Empty transition"
        return (
            f'  {nodes[head]} -> {nodes[tail]} '
            f'[label="&lambda;",style=dashed,'
            f'tooltip="{empty_tooltip}",edgetooltip="{empty_tooltip}",labeltooltip="{empty_tooltip}"]'
        )


def _style_nodes(nodes, node_info, source_node, sink_node, dominator_nodes):
    """
    Generate DOT notation for node styling (source, sink, dominators, regular nodes).
    """
    dot_lines = []

    # Build set of all styled nodes
    styled_nodes = set()
    if sink_node:
        styled_nodes.add(sink_node)
    if source_node:
        styled_nodes.add(source_node)
    styled_nodes.update(dominator_nodes)

    # Style sink node
    if sink_node:
        dot_lines.append(
            f'{nodes[sink_node]}'
            f'[fillcolor=aliceblue,style=filled,peripheries=2,penwidth=2,'
            f'tooltip="Sink node\\n{node_info[sink_node]}"]'
        )

    # Style source node
    if source_node:
        dot_lines.extend([
            f'{nodes[source_node]}'
            f'[fillcolor=aliceblue,style=filled,penwidth=2,'
            f'tooltip="Source node\\n{node_info[source_node]}"]',
            'i[label="",shape=point,width=.1,tooltip="Start"]',
            f'i->{nodes[source_node]}',
        ])

    # Style dominator nodes (excluding source/sink to avoid overwriting)
    for node in dominator_nodes:
        if node not in [sink_node, source_node]:
            dot_lines.append(
                f'{nodes[node]}'
                f'[fillcolor=aliceblue,style=filled,penwidth=2,'
                f'tooltip="Dominator node\\n{node_info[node]}"]'
            )

    # Add tooltips for regular nodes (those without special styling)
    for node, node_id in nodes.items():
        if node not in styled_nodes:
            dot_lines.append(f'{node_id}[tooltip="{node_info[node]}"]')

    return dot_lines


def graph_to_dot(graph, reference, selector=None, dominators=True, edges_limit=100):
    """
    Convert a graph to DOT format for visualization purposes.

    It includes styled nodes (source, sink, dominators) and edges with variants.
    Tooltips are added to both nodes and edges for additional information.
    The source node is always labeled as s0, and the sink node gets the highest number.
    """
    dot_lines = [
        "digraph {",
        "rankdir=LR",
        "edge[fontname=monospace]",
        f'node[fixedsize=true,fontname=serif,shape=circle,width=1]'
    ]

    summary = _collect_nodes(graph, selector, edges_limit)
    if summary[4]:
        return f"// Graph too complex. Stopped at {summary[3]} edges and {len(summary[0])} nodes."

    all_nodes, head_nodes, tail_nodes, _, _ = summary

    source_node = next(iter(head_nodes - tail_nodes), None)
    sink_node = next(iter(tail_nodes - head_nodes), None)

    nodes, node_info = _assign_node_ids(all_nodes, source_node, sink_node)

    for edge in graph.edges():
        head, tail, variant, count = edge["head"], edge["tail"], edge["variant"], edge["count"]
        if selector and selector[-1] is True:
            head, tail = tail, head
        edge_dot = _create_edge_dot(head, tail, variant, count, reference, selector, nodes)
        dot_lines.append(edge_dot)

    dominator_nodes = get_dominators(graph) if dominators else set()
    node_styling = _style_nodes(nodes, node_info, source_node, sink_node, dominator_nodes)
    dot_lines.extend(node_styling)
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



def to_hgvs(variant, reference=None, selector=None, only_substitutions=True):
    """
    Adapted from the algebra.
    """
    def _seq(sequence):
        if selector is not None and selector[-1]:
            return reverse_complement(sequence)
        return sequence

    def _loc(start, end=None):
        if selector is not None:
            if end is not None:
                if selector[-1]:
                    return f"{g_to_cn(end, selector)}_{g_to_cn(start, selector)}"
                return f"{g_to_cn(start, selector)}_{g_to_cn(end, selector)}"
            return f"{g_to_cn(start, selector)}"
        if end:
            return f"{start}_{end}"
        return start

    if variant.end - variant.start == 0:
        if not variant.sequence:
            return "="
        return f"{_loc(variant.start, variant.start + 1)}ins{_seq(variant.sequence)}"

    deleted = ""
    substitution = ""
    if reference is not None:
        if not only_substitutions:
            deleted = reference[variant.start:variant.end]
        substitution = reference[variant.start:variant.end]

    if variant.end - variant.start == 1:
        if not variant.sequence:
            return f"{_loc(variant.start + 1)}del{deleted}"
        if len(variant.sequence) == 1:
            return f"{_loc(variant.start + 1)}{substitution}>{_seq(variant.sequence)}"
        return f"{_loc(variant.start + 1)}del{deleted}ins{_seq(variant.sequence)}"

    if not variant.sequence:
        return f"{_loc(variant.start + 1, variant.end)}del{deleted}"

    return f"{_loc(variant.start + 1, variant.end)}del{deleted}ins{_seq(variant.sequence)}"


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
