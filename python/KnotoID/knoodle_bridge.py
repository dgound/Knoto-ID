"""Interoperability bridge between KnotoID and Knoodle (pyknoodle).

Both tools represent a planar diagram with a 0-based PD code using the
KnotTheory column order [under_in, over, under_out, over_other]. The only
difference is the crossing sign: Knoodle carries it as an explicit 5th
"handedness" column (+1 = right, -1 = left), whereas KnotoID encodes it
implicitly in the order of the four arc labels. The bridge therefore just
drops (import) or omits (export) that column.

The ``knoodle`` package is only needed for :func:`to_knoodle`; it is imported
lazily so that :func:`from_knoodle` and the rest of KnotoID work without it.

Examples
--------
    from KnotoID.knoodle_bridge import from_knoodle, to_knoodle
    import knoodle

    ka = knoodle.KnotAnalyzer(coordinates)     # a Knoodle diagram
    diagram = from_knoodle(ka)                  # -> knotoID_cpp.PlanarDiagram
    ka2 = to_knoodle(diagram)                   # KnotoID -> Knoodle
"""

import knotoID_cpp as _kn


def from_knoodle(source, flag_planar=False, flag_debug=False):
    """Build a KnotoID ``PlanarDiagram`` from Knoodle PD data.

    :param source: either a Knoodle ``KnotAnalyzer`` (anything exposing
        ``get_pd_code_matrix()``) or a plain PD matrix -- a list of
        ``[a0, a1, a2, a3]`` or ``[a0, a1, a2, a3, handedness]`` rows.
    :param flag_planar: build a planar (plane) diagram rather than a spherical
        one. Note: planar open knotoids also need exterior arcs, which a bare
        Knoodle PD code does not carry.
    :param flag_debug: enable the underlying C++ debug output.
    :returns: a ``knotoID_cpp.PlanarDiagram``.
    """
    matrix = source.get_pd_code_matrix() if hasattr(source, "get_pd_code_matrix") else source
    # keep the four arc labels; drop Knoodle's handedness column if present
    crossings = [list(row)[:4] for row in matrix]
    diagram = _kn.PlanarDiagram(flag_planar, flag_debug)
    diagram.load_from_pd_code(crossings)
    return diagram


def to_knoodle(diagram, simplify=True, simplify_level=5):
    """Send a KnotoID ``PlanarDiagram`` to Knoodle as a ``KnotAnalyzer``.

    :param diagram: a ``knotoID_cpp.PlanarDiagram``.
    :param simplify: passed through to ``KnotAnalyzer.from_pd_code``.
    :param simplify_level: passed through to ``KnotAnalyzer.from_pd_code``.
    :returns: a ``knoodle.KnotAnalyzer`` built from the diagram's PD code.
    :raises ImportError: if the ``knoodle`` package is not installed.
    """
    import knoodle
    return knoodle.KnotAnalyzer.from_pd_code(
        diagram.to_pd_code(), simplify=simplify, simplify_level=simplify_level)
