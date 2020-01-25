"""
The :mod:`scadnano` Python module is a library for describing synthetic DNA nanostructures
(e.g., DNA origami).
Installation instructions are at the
`GitHub repository <https://github.com/UC-Davis-molecular-computing/scadnano-python-package>`_.

This module is used to write Python scripts outputting ``*.dna`` files readable
by `scadnano <https://web.cs.ucdavis.edu/~doty/scadnano/>`_, a web application useful for displaying
and manually editing these structures.
The purpose of this module is to help automate some of the task of creating DNA designs,
as well as making large-scale changes to them that are easier to describe programmatically than
to do by hand in scadnano.

This library uses typing hints from the Python typing library.
(https://docs.python.org/3/library/typing.html)
Each function and method indicate intended types of the parameters.
However, due to Python's design, these types are not enforced at runtime.
It is suggested to use a static analysis tool such as that provided by an IDE such as PyCharm
(https://www.jetbrains.com/pycharm/)
to see warnings when the typing rules are violated. 
Such warnings probably indicate an erroneous usage.

Most of the classes in this module are Python dataclasses
(https://docs.python.org/3/library/dataclasses.html)
whose fields show up in the documentation.
Their types are listed in parentheses after the name of the class;
for example :any:`Color` has ``int`` fields :py:data:`Color.r`, :py:data:`Color.g`, :py:data:`Color.b`.
In general it is safe to read these fields directly, but not to write to them directly.
Setter methods (named ``set_<fieldname>``) are provided for fields where it makes sense to set it to another
value than it had originally.
However, due to Python naming conventions for dataclass fields and property setters,
it is not straightforward to enforce that the fields cannot be written, 
so the user must take care not to set them.
"""

# needed to use forward annotations: https://docs.python.org/3/whatsnew/3.7.html#whatsnew37-pep563
from __future__ import annotations

from abc import abstractmethod, ABC
import json
import enum
import itertools
import re
from dataclasses import dataclass, field, InitVar
from typing import Tuple, List, Dict, Union, Optional
from collections import defaultdict, OrderedDict
import sys
import os.path
import xlwt
from docutils.nodes import subscript


# from scadnano.json_utils import JSONSerializable, json_encode, NoIndent


# TODO: add Boolean field Strand.circular

# TODO: make explicit rules about when strands can be added and sequences assigned.
#  For instance, if we add a strand to overlap one that already has a DNA sequence sequence assigned,
#  should the complement be automatically assigned?

# TODO: add support for writing 3D positions (in addition to 2D svg_positions)

# TODO: add support for writing files uploadable to other synthesis company websites besides IDT

# TODO: see if :param the_paramter: and :return: can be used with Sphinx


def _pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


##############################################################################
# JSON serialization
# There are external libraries to handle JSON
# in Python, but I want this to be a simple, single-file library, so we just
# implement what we need below.


class _JSONSerializable(ABC):

    @abstractmethod
    def to_json_serializable(self, suppress_indent=True):
        raise NotImplementedError()


def _json_encode(obj: _JSONSerializable, suppress_indent: bool = True) -> str:
    encoder = _SuppressableIndentEncoder if suppress_indent else json.JSONEncoder
    serializable = obj.to_json_serializable(suppress_indent=suppress_indent)
    return json.dumps(serializable, cls=encoder, indent=2)


class _NoIndent:
    """ Value wrapper. """

    def __init__(self, value):
        self.value = value


class _SuppressableIndentEncoder(json.JSONEncoder):
    def __init__(self, *args, **kwargs):
        self.unique_id = 0
        super(_SuppressableIndentEncoder, self).__init__(*args, **kwargs)
        self.kwargs = dict(kwargs)
        del self.kwargs['indent']
        self._replacement_map = {}

    def default(self, obj):
        if isinstance(obj, _NoIndent):
            # key = uuid.uuid1().hex # this caused problems with Brython.
            key = self.unique_id
            self.unique_id += 1
            self._replacement_map[key] = json.dumps(obj.value, **self.kwargs)
            return "@@%s@@" % (key,)
        else:
            return super().default(obj)

    def encode(self, obj):
        result = super().encode(obj)
        for k, v in self._replacement_map.items():
            result = result.replace('"@@%s@@"' % (k,), v)
        return result


#
# END JSON serialization
##############################################################################


##############################################################################
# Colors
# As with JSON serialization, there are external libraries to handle colors
# in Python, but I want this to be a simple, single-file library, so we just
# implement what we need below.

@dataclass
class Color(_JSONSerializable):
    r: int = None
    """Red component: 0-255.
    
    Optional if :py:data:`Color.hex` is given."""

    g: int = None
    """Green component: 0-255.
    
    Optional if :py:data:`Color.hex` is given."""

    b: int = None
    """Blue component: 0-255.
    
    Optional if :py:data:`Color.hex` is given."""

    hex: InitVar[str] = None
    """Hex color preceded by # sign, e.g., "#ff0000" is red.
    
    Optional if :py:data:`Color.r`, :py:data:`Color.g`, :py:data:`Color.b` are all given."""

    def __post_init__(self, hex):
        if hex is None:
            assert (self.r is not None and self.g is not None and self.b is not None)
        else:
            assert (self.r is None and self.g is None and self.b is None)
            hex = hex.lstrip('#')
            self.r = int(hex[0:2], 16)
            self.g = int(hex[2:4], 16)
            self.b = int(hex[4:6], 16)

    def to_json_serializable(self, suppress_indent=True):
        # Return object representing this Color that is JSON serializable.
        # return NoIndent(self.__dict__) if suppress_indent else self.__dict__
        return f'#{self.r:02x}{self.g:02x}{self.b:02x}'

    def to_cadnano_v2_int_hex(self):
        return int(f'{self.r:02x}{self.g:02x}{self.b:02x}', 16)


class ColorCycler:
    """
    Calling ``next(color_cycler)`` on a ColorCycler named ``color_cycler``
    returns a the next :any:`Color` from a fixed size list,
    cycling after reaching the end of the list.

    To choose new colors, set ``color_cycler.colors`` to a new list of :any:`Color`'s.
    """

    # These are copied from cadnano:
    # https://github.com/sdouglas/cadnano2/blob/master/views/styles.py#L97
    _colors = [Color(50, 184, 108),
               Color(204, 0, 0),
               Color(247, 67, 8),
               Color(247, 147, 30),
               Color(170, 170, 0),
               Color(87, 187, 0),
               Color(0, 114, 0),
               Color(3, 182, 162),
               # Color(23, 0, 222), # don't like this because it looks too much like scaffold
               Color(50, 0, 150),  # this one is better contrast with scaffold
               Color(184, 5, 108),
               Color(51, 51, 51),
               Color(115, 0, 222),
               Color(136, 136, 136)]
    """List of colors to cycle through."""

    def __init__(self):
        self._current_color_idx = 0
        # random order
        order = [3, 11, 0, 12, 8, 1, 10, 6, 5, 9, 4, 7, 2]
        colors_shuffled = [None] * len(self._colors)
        for i, color in zip(order, self._colors):
            colors_shuffled[i] = color
        self._colors = colors_shuffled

    def __next__(self):
        color = self._colors[self._current_color_idx]
        self._current_color_idx = (self._current_color_idx + 1) % len(self._colors)
        return color

    @property
    def colors(self):
        """The colors that are cycled through when calling ``next()`` on some :any:`ColorCycler`."""
        return list(self._colors)

    @colors.setter
    def colors(self, newcolors):
        self._colors = newcolors
        self._current_color_idx = 0


default_scaffold_color = Color(0, 102, 204)
"""Default color for scaffold strand(s)."""

default_strand_color = Color(0, 0, 0)
"""Default color for non-scaffold strand(s)."""

color_cycler = ColorCycler()


#
# END Colors
##############################################################################


@enum.unique
class Grid(str, enum.Enum):
    """
    Represents default patterns for laying out helices in the side view.
    Each :any:`Grid` except :py:data:`Grid.none` has an interpretation of a "grid position",
    which is a 2D integer coordinate (`h`, `v`).
    (scadnano also allows a 3rd coordinate (`h`, `v`, `b`) specifying a "base offset" at which to position
    the start of the :any:`Helix`, which is not relevant for the side view but will eventually be
    supported to adjust the main view.)
    """

    square = "square"
    """
    Square lattice. 
    Increasing `h` moves right and increasing `v` moves down. 
    (i.e., "computer screen coordinates" rather than Cartesian coordinates where positive `y` is up.)
    """

    hex = "hex"
    """
    Hexagonal lattice. Uses the *"odd-r horizontal layout"* coordinate system described here: 
    https://www.redblobgames.com/grids/hexagons/. 
    Incrementing `v` moves down and to the right if `h` is even, 
    and moves down and to the left if `h` is odd.
    """

    honeycomb = "honeycomb"
    """
    Honeycomb lattice. This consists of all the hex lattice positions except where 
    honeycomb lattice disallows grid positions (`h`, `v`) with 
    `v` even and `h` a multiple of 3 or
    `v` odd and `h` = 1 + a multiple of 3.
    """

    none = "none"
    """No fixed grid."""


# convenience names for users
reverse = False
forward = True
square = Grid.square
hexagonal = Grid.hex  # should not use identifier "hex" because that's a Python built-in function
honeycomb = Grid.honeycomb

##########################################################################
# constants

current_version: str = "0.1.0"
initial_version: str = "0.1.0"

default_idt_scale = "25nm"
default_idt_purification = "STD"


def default_major_tick_distance(grid: Grid) -> int:
    return 7 if grid in (Grid.hex, Grid.honeycomb) else 8


default_grid: Grid = Grid.square

default_helix_rotation: float = -90.0
default_helix_rotation_anchor: int = 0

base_width_svg: float = 10.0
"""Width of a single base in the SVG main view of scadnano."""

base_height_svg: float = 10.0
"""Height of a single base in the SVG main view of scadnano."""

distance_between_helices_svg: float = (base_width_svg * 2.5 / 0.34)
"""Distance between tops of two consecutive helices (using default positioning rules).

This is set to (:const:`base_width_svg` * 2.5/0.34) based on the following calculation,
to attempt to make the DNA structure appear to scale in 2D drawings:
The width of one base pair of double-stranded DNA bp is 0.34 nm. In a DNA origami, 
AFM images let us estimate that the average distance between adjacent double helices is 2.5 nm.
(A DNA double-helix is only 2 nm wide, but the helices electrostatically repel each other so the spacing
in a DNA origami or an other DNA nanostructure with many parallel DNA helices---e.g., single-stranded tile
lattices---is larger than 2 nm.)
Thus the distance between the helices is 2.5/0.34 ~ 7.5 times the width of a single DNA base.
"""

DNA_base_wildcard: str = '?'
"""Symbol to insert when a DNA sequence has been assigned to a strand through complementarity, but
some regions of the strand are not bound to the strand that was just assigned. Also used in case the
DNA sequence assigned to a strand is too short; the sequence is padded with :any:`DNA_base_wildcard` to 
make its length the same as the length of the strand."""

m13_sequence = \
    "TTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACC" \
    "CCAAAAAACTTGATTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGG" \
    "ACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGAACCACCATCAAACAGGATT" \
    "TTCGCCTGCTGGGGCAAACCAGCGTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCCGTCTCACTGGTGAAAAG" \
    "AAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGG" \
    "CAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGC" \
    "GGATAACAATTTCACACAGGAAACAGCTATGACCATGATTACGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTG" \
    "GCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATA" \
    "GCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAG" \
    "CTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTAT" \
    "CCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGA" \
    "CGCGAATTATTTTTGATGGCGTTCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAATGCGAATTTTAACAAAATATTAACGTTTACAATTTAA" \
    "ATATTTGCTTATACAATCTTCCTGTTTTTGGGGCTTTTCTGATTATCAACCGGGGTACATATGATTGACATGCTAGTTTTACGATTACCGTTCATCGATT" \
    "CTCTTGTTTGCTCCAGACTCTCAGGCAATGACCTGATAGCCTTTGTAGATCTCTCAAAAATAGCTACCCTCTCCGGCATTAATTTATCAGCTAGAACGGT" \
    "TGAATATCATATTGATGGTGATTTGACTGTCTCCGGCCTTTCTCACCCTTTTGAATCTTTACCTACACATTACTCAGGCATTGCATTTAAAATATATGAG" \
    "GGTTCTAAAAATTTTTATCCTTGCGTTGAAATAAAGGCTTCTCCCGCAAAAGTATTACAGGGTCATAATGTTTTTGGTACAACCGATTTAGCTTTATGCT" \
    "CTGAGGCTTTATTGCTTAATTTTGCTAATTCTTTGCCTTGCCTGTATGATTTATTGGATGTTAATGCTACTACTATTAGTAGAATTGATGCCACCTTTTC" \
    "AGCTCGCGCCCCAAATGAAAATATAGCTAAACAGGTTATTGACCATTTGCGAAATGTATCTAATGGTCAAACTAAATCTACTCGTTCGCAGAATTGGGAA" \
    "TCAACTGTTATATGGAATGAAACTTCCAGACACCGTACTTTAGTTGCATATTTAAAACATGTTGAGCTACAGCATTATATTCAGCAATTAAGCTCTAAGC" \
    "CATCCGCAAAAATGACCTCTTATCAAAAGGAGCAATTAAAGGTACTCTCTAATCCTGACCTGTTGGAGTTTGCTTCCGGTCTGGTTCGCTTTGAAGCTCG" \
    "AATTAAAACGCGATATTTGAAGTCTTTCGGGCTTCCTCTTAATCTTTTTGATGCAATCCGCTTTGCTTCTGACTATAATAGTCAGGGTAAAGACCTGATT" \
    "TTTGATTTATGGTCATTCTCGTTTTCTGAACTGTTTAAAGCATTTGAGGGGGATTCAATGAATATTTATGACGATTCCGCAGTATTGGACGCTATCCAGT" \
    "CTAAACATTTTACTATTACCCCCTCTGGCAAAACTTCTTTTGCAAAAGCCTCTCGCTATTTTGGTTTTTATCGTCGTCTGGTAAACGAGGGTTATGATAG" \
    "TGTTGCTCTTACTATGCCTCGTAATTCCTTTTGGCGTTATGTATCTGCATTAGTTGAATGTGGTATTCCTAAATCTCAACTGATGAATCTTTCTACCTGT" \
    "AATAATGTTGTTCCGTTAGTTCGTTTTATTAACGTAGATTTTTCTTCCCAACGTCCTGACTGGTATAATGAGCCAGTTCTTAAAATCGCATAAGGTAATT" \
    "CACAATGATTAAAGTTGAAATTAAACCATCTCAAGCCCAATTTACTACTCGTTCTGGTGTTTCTCGTCAGGGCAAGCCTTATTCACTGAATGAGCAGCTT" \
    "TGTTACGTTGATTTGGGTAATGAATATCCGGTTCTTGTCAAGATTACTCTTGATGAAGGTCAGCCAGCCTATGCGCCTGGTCTGTACACCGTTCATCTGT" \
    "CCTCTTTCAAAGTTGGTCAGTTCGGTTCCCTTATGATTGACCGTCTGCGCCTCGTTCCGGCTAAGTAACATGGAGCAGGTCGCGGATTTCGACACAATTT" \
    "ATCAGGCGATGATACAAATCTCCGTTGTACTTTGTTTCGCGCTTGGTATAATCGCTGGGGGTCAAAGATGAGTGTTTTAGTGTATTCTTTTGCCTCTTTC" \
    "GTTTTAGGTTGGTGCCTTCGTAGTGGCATTACGTATTTTACCCGTTTAATGGAAACTTCCTCATGAAAAAGTCTTTAGTCCTCAAAGCCTCTGTAGCCGT" \
    "TGCTACCCTCGTTCCGATGCTGTCTTTCGCTGCTGAGGGTGACGATCCCGCAAAAGCGGCCTTTAACTCCCTGCAAGCCTCAGCGACCGAATATATCGGT" \
    "TATGCGTGGGCGATGGTTGTTGTCATTGTCGGCGCAACTATCGGTATCAAGCTGTTTAAGAAATTCACCTCGAAAGCAAGCTGATAAACCGATACAATTA" \
    "AAGGCTCCTTTTGGAGCCTTTTTTTTGGAGATTTTCAACGTGAAAAAATTATTATTCGCAATTCCTTTAGTTGTTCCTTTCTATTCTCACTCCGCTGAAA" \
    "CTGTTGAAAGTTGTTTAGCAAAATCCCATACAGAAAATTCATTTACTAACGTCTGGAAAGACGACAAAACTTTAGATCGTTACGCTAACTATGAGGGCTG" \
    "TCTGTGGAATGCTACAGGCGTTGTAGTTTGTACTGGTGACGAAACTCAGTGTTACGGTACATGGGTTCCTATTGGGCTTGCTATCCCTGAAAATGAGGGT" \
    "GGTGGCTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTTCTGAGGGTGGCGGTACTAAACCTCCTGAGTACGGTGATACACCTATTCCGGGCTATACTTATA" \
    "TCAACCCTCTCGACGGCACTTATCCGCCTGGTACTGAGCAAAACCCCGCTAATCCTAATCCTTCTCTTGAGGAGTCTCAGCCTCTTAATACTTTCATGTT" \
    "TCAGAATAATAGGTTCCGAAATAGGCAGGGGGCATTAACTGTTTATACGGGCACTGTTACTCAAGGCACTGACCCCGTTAAAACTTATTACCAGTACACT" \
    "CCTGTATCATCAAAAGCCATGTATGACGCTTACTGGAACGGTAAATTCAGAGACTGCGCTTTCCATTCTGGCTTTAATGAGGATTTATTTGTTTGTGAAT" \
    "ATCAAGGCCAATCGTCTGACCTGCCTCAACCTCCTGTCAATGCTGGCGGCGGCTCTGGTGGTGGTTCTGGTGGCGGCTCTGAGGGTGGTGGCTCTGAGGG" \
    "TGGCGGTTCTGAGGGTGGCGGCTCTGAGGGAGGCGGTTCCGGTGGTGGCTCTGGTTCCGGTGATTTTGATTATGAAAAGATGGCAAACGCTAATAAGGGG" \
    "GCTATGACCGAAAATGCCGATGAAAACGCGCTACAGTCTGACGCTAAAGGCAAACTTGATTCTGTCGCTACTGATTACGGTGCTGCTATCGATGGTTTCA" \
    "TTGGTGACGTTTCCGGCCTTGCTAATGGTAATGGTGCTACTGGTGATTTTGCTGGCTCTAATTCCCAAATGGCTCAAGTCGGTGACGGTGATAATTCACC" \
    "TTTAATGAATAATTTCCGTCAATATTTACCTTCCCTCCCTCAATCGGTTGAATGTCGCCCTTTTGTCTTTGGCGCTGGTAAACCATATGAATTTTCTATT" \
    "GATTGTGACAAAATAAACTTATTCCGTGGTGTCTTTGCGTTTCTTTTATATGTTGCCACCTTTATGTATGTATTTTCTACGTTTGCTAACATACTGCGTA" \
    "ATAAGGAGTCTTAATCATGCCAGTTCTTTTGGGTATTCCGTTATTATTGCGTTTCCTCGGTTTCCTTCTGGTAACTTTGTTCGGCTATCTGCTTACTTTT" \
    "CTTAAAAAGGGCTTCGGTAAGATAGCTATTGCTATTTCATTGTTTCTTGCTCTTATTATTGGGCTTAACTCAATTCTTGTGGGTTATCTCTCTGATATTA" \
    "GCGCTCAATTACCCTCTGACTTTGTTCAGGGTGTTCAGTTAATTCTCCCGTCTAATGCGCTTCCCTGTTTTTATGTTATTCTCTCTGTAAAGGCTGCTAT" \
    "TTTCATTTTTGACGTTAAACAAAAAATCGTTTCTTATTTGGATTGGGATAAATAATATGGCTGTTTATTTTGTAACTGGCAAATTAGGCTCTGGAAAGAC" \
    "GCTCGTTAGCGTTGGTAAGATTCAGGATAAAATTGTAGCTGGGTGCAAAATAGCAACTAATCTTGATTTAAGGCTTCAAAACCTCCCGCAAGTCGGGAGG" \
    "TTCGCTAAAACGCCTCGCGTTCTTAGAATACCGGATAAGCCTTCTATATCTGATTTGCTTGCTATTGGGCGCGGTAATGATTCCTACGATGAAAATAAAA" \
    "ACGGCTTGCTTGTTCTCGATGAGTGCGGTACTTGGTTTAATACCCGTTCTTGGAATGATAAGGAAAGACAGCCGATTATTGATTGGTTTCTACATGCTCG" \
    "TAAATTAGGATGGGATATTATTTTTCTTGTTCAGGACTTATCTATTGTTGATAAACAGGCGCGTTCTGCATTAGCTGAACATGTTGTTTATTGTCGTCGT" \
    "CTGGACAGAATTACTTTACCTTTTGTCGGTACTTTATATTCTCTTATTACTGGCTCGAAAATGCCTCTGCCTAAATTACATGTTGGCGTTGTTAAATATG" \
    "GCGATTCTCAATTAAGCCCTACTGTTGAGCGTTGGCTTTATACTGGTAAGAATTTGTATAACGCATATGATACTAAACAGGCTTTTTCTAGTAATTATGA" \
    "TTCCGGTGTTTATTCTTATTTAACGCCTTATTTATCACACGGTCGGTATTTCAAACCATTAAATTTAGGTCAGAAGATGAAATTAACTAAAATATATTTG" \
    "AAAAAGTTTTCTCGCGTTCTTTGTCTTGCGATTGGATTTGCATCAGCATTTACATATAGTTATATAACCCAACCTAAGCCGGAGGTTAAAAAGGTAGTCT" \
    "CTCAGACCTATGATTTTGATAAATTCACTATTGACTCTTCTCAGCGTCTTAATCTAAGCTATCGCTATGTTTTCAAGGATTCTAAGGGAAAATTAATTAA" \
    "TAGCGACGATTTACAGAAGCAAGGTTATTCACTCACATATATTGATTTATGTACTGTTTCCATTAAAAAAGGTAATTCAAATGAAATTGTTAAATGTAAT" \
    "TAATTTTGTTTTCTTGATGTTTGTTTCATCATCTTCTTTTGCTCAGGTAATTGAAATGAATAATTCGCCTCTGCGCGATTTTGTAACTTGGTATTCAAAG" \
    "CAATCAGGCGAATCCGTTATTGTTTCTCCCGATGTAAAAGGTACTGTTACTGTATATTCATCTGACGTTAAACCTGAAAATCTACGCAATTTCTTTATTT" \
    "CTGTTTTACGTGCAAATAATTTTGATATGGTAGGTTCTAACCCTTCCATTATTCAGAAGTATAATCCAAACAATCAGGATTATATTGATGAATTGCCATC" \
    "ATCTGATAATCAGGAATATGATGATAATTCCGCTCCTTCTGGTGGTTTCTTTGTTCCGCAAAATGATAATGTTACTCAAACTTTTAAAATTAATAACGTT" \
    "CGGGCAAAGGATTTAATACGAGTTGTCGAATTGTTTGTAAAGTCTAATACTTCTAAATCCTCAAATGTATTATCTATTGACGGCTCTAATCTATTAGTTG" \
    "TTAGTGCTCCTAAAGATATTTTAGATAACCTTCCTCAATTCCTTTCAACTGTTGATTTGCCAACTGACCAGATATTGATTGAGGGTTTGATATTTGAGGT" \
    "TCAGCAAGGTGATGCTTTAGATTTTTCATTTGCTGCTGGCTCTCAGCGTGGCACTGTTGCAGGCGGTGTTAATACTGACCGCCTCACCTCTGTTTTATCT" \
    "TCTGCTGGTGGTTCGTTCGGTATTTTTAATGGCGATGTTTTAGGGCTATCAGTTCGCGCATTAAAGACTAATAGCCATTCAAAAATATTGTCTGTGCCAC" \
    "GTATTCTTACGCTTTCAGGTCAGAAGGGTTCTATCTCTGTTGGCCAGAATGTCCCTTTTATTACTGGTCGTGTGACTGGTGAATCTGCCAATGTAAATAA" \
    "TCCATTTCAGACGATTGAGCGTCAAAATGTAGGTATTTCCATGAGCGTTTTTCCTGTTGCAATGGCTGGCGGTAATATTGTTCTGGATATTACCAGCAAG" \
    "GCCGATAGTTTGAGTTCTTCTACTCAGGCAAGTGATGTTATTACTAATCAAAGAAGTATTGCTACAACGGTTAATTTGCGTGATGGACAGACTCTTTTAC" \
    "TCGGTGGCCTCACTGATTATAAAAACACTTCTCAGGATTCTGGCGTACCGTTCCTGTCTAAAATCCCTTTAATCGGCCTCCTGTTTAGCTCCCGCTCTGA" \
    "TTCTAACGAGGAAAGCACGTTATACGTGCTCGTCAAAGCAACCATAGTACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAG" \
    "CGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTC"
"""
The M13mp18 DNA sequence (commonly called simply M13), starting from cyclic rotation 5588, as defined in
`GenBank <https://www.neb.com/~/media/NebUs/Page%20Images/Tools%20and%20Resources/Interactive%20Tools/DNA%20Sequences%20and%20Maps/Text%20Documents/m13mp18gbk.txt>`_.
This is the "standard" variant of consisting of 7249 bases, sold by companies such as  
`New England Biolabs <https://www.neb.com/~/media/nebus/page%20images/tools%20and%20resources/interactive%20tools/dna%20sequences%20and%20maps/m13mp18_map.pdf>`_
and
`Tilibit <https://cdn.shopify.com/s/files/1/1299/5863/files/Product_Sheet_single-stranded_scaffold_DNA_type_7249_M1-10.pdf?14656642867652657391>`_.

The actual M13 DNA strand itself is circular, 
so assigning this sequence to the scaffold :any:`Strand` in a :any:`DNADesign`
means that the "5' end" of the scaffold :any:`Strand` 
(which is a fiction since the actual circular DNA strand has no endpoint) 
will have the sequence starting at position 5588 starting at the displayed 5' in scadnano,
assigned until the displayed 3' end. 
Assuming the displayed scaffold :any:`Strand` has length :math:`n < 7249`, then a loopout of length 
:math:`7249 - n` consisting of the undisplayed bases will be present in the actual DNA structure.
For a more detailed discussion of why this particular rotation of M13 is chosen,
see 
`Supplementary Note S8 <http://www.dna.caltech.edu/Papers/DNAorigami-supp1.linux.pdf>`_ 
in
[`Folding DNA to create nanoscale shapes and patterns. Paul W. K. Rothemund, Nature 440:297-302 (2006) <http://www.nature.com/nature/journal/v440/n7082/abs/nature04586.html>`_].
"""  # noqa (suppress PEP warning)

##################
# keys

# DNADesign keys
version_key = 'version'
grid_key = 'grid'
major_tick_distance_key = 'major_tick_distance'
major_ticks_key = 'major_ticks'
rotation_key = 'rotation'
rotation_anchor_key = 'rotation_anchor'
helices_key = 'helices'
strands_key = 'strands'
scaffold_key = 'scaffold'
helices_view_order_key = "helices_view_order"
is_origami_key = 'is_origami'

# Helix keys
idx_key = 'idx'
max_offset_key = 'max_offset'
min_offset_key = 'min_offset'
grid_position_key = 'grid_position'
svg_position_key = 'svg_position'
position3d_key = 'position'

# Strand keys
color_key = 'color'
dna_sequence_key = 'dna_sequence'
substrands_key = 'substrands'
idt_key = 'idt'
is_scaffold_key = 'is_scaffold'

# Substrand keys
helix_idx_key = 'helix'
forward_key = 'forward'
start_key = 'start'
end_key = 'end'
deletions_key = 'deletions'
insertions_key = 'insertions'

# Loopout keys
loopout_key = 'loopout'


# end keys
##################

# end constants
##########################################################################

@dataclass
class Position3D(_JSONSerializable):
    """
    Position (x,y,z) and orientation (pitch,roll,yaw) in 3D space.
    """

    x: float
    y: float
    z: float
    pitch: float
    roll: float
    yaw: float

    def to_json_serializable(self, suppress_indent=True):
        dct = self.__dict__
        # return NoIndent(dct) if suppress_indent else dct
        return dct


def in_browser() -> bool:
    """Test if this code is running in the browser.

    Checks for existence of package "pyodide" used in pyodide. If present it is assumed the code is
    running in the browser."""
    try:
        import pyodide
        return True
    except ModuleNotFoundError:
        return False


# TODO: add rotation field to Helix
#   doc: https://docs.google.com/document/d/1OysNEI1RIwzJ6IqbfTgLR3fYsmEUqIiuHJKxO7CeuaQ/edit#

@dataclass
class Helix(_JSONSerializable):
    """
    Represents a "helix" where :any:`Substrand`'s could go. Technically a :any:`Helix` can contain no
    :any:`Substrand`'s. More commonly, some partial regions of it may have only 1 or 0 :any:`Substrand`'s.
    So it is best thought of as a "potential" double-helix.

    It has a 1-dimensional integer coordinate system given by "offsets", integers between
    :py:data:`Helix.min_offset` (inclusive) and :py:data:`Helix.max_offset` (exclusive).
    At any valid offset for this :any:`Helix`, at most two :any:`Substrand`'s may share that offset
    on this :any:`Helix`, and if there are exactly two, then one must have
    :py:data:`Substrand.forward` = ``true`` and the other must have
    :py:data:`Substrand.forward` = ``false``.

    Once part of a :any:`DNADesign`, a :any:`Helix` has an index (accessible  via :py:meth:`Helix.idx`
    once the :any:`DNADesign` is created)
    representing its order in the list of all :any:`Helix`'s. This index is how a :any:`Substrand` is
    associated to the :any:`Helix` via the integer index :any:`Substrand.helix`.
    """

    max_offset: int = None
    """Maximum offset (exclusive) of :any:`Substrand` that can be drawn on this :any:`Helix`. 
    If unspecified, it is calculated when the :any:`DNADesign` is instantiated as 
    the largest :any:`Substrand.end` offset of any :any:`Substrand` in the design.
    """

    min_offset: int = 0
    """Minimum offset (inclusive) of :any:`Substrand` that can be drawn on this :any:`Helix`. 
    If unspecified, it is set to 0.
    """

    major_tick_distance: int = -1
    """If positive, overrides :any:`DNADesign.major_tick_distance`."""

    major_ticks: List[int] = None
    """If not ``None``, overrides :any:`DNADesign.major_tick_distance` and :any:`Helix.major_tick_distance`
    to specify a list of offsets at which to put major ticks."""

    grid_position: Tuple[int, int, int] = None
    """`(h,v,b)` position of this helix in the side view grid,
    if :const:`Grid.square`, :const:`Grid.hex` , or :const:`Grid.honeycomb` is used
    in the :any:`DNADesign` containing this helix.
    `h` and `v` are in units of "helices": incrementing `h` moves right one helix in the grid
    and incrementing `v` moves down one helix in the grid. 
    In the case of the hexagonal or honeycomb lattice, 
    The convention is that incrementing `v` moves down and to the right if h is even, 
    and moves down and to the left if `h` is odd.
    This is the "odd-r horizontal layout" coordinate system here: 
    https://www.redblobgames.com/grids/hexagons/)
    `b` goes in and out of the screen in the side view, and it is in units of "bases".
    Incrementing `b` moves the whole helix one base into the screen.
    In the main view, a helix with `b` = 1 would have its base offset 0 line up with base offset 1
    of a helix with `b` = 0.
    However, the default y svg_position for helices does not otherwise depend on grid_position.
    The default is to list the y-coordinates in order by helix idx.
    
    Default is `h` = 0, `v` = index of :any:`Helix` in :py:data:`DNADesign.helices`, `b` = 0."""

    svg_position: Tuple[float, float] = None
    """`(x,y)` SVG coordinates of base offset 0 of this Helix in the main view. 
    
    If `grid_position` and `position` are both omitted, then the default is 
    `x` = 0, `y` = [index of helix] * :any:`scadnano.distance_between_helices_svg`.
    
    If `grid_position = (h,v,b)` is specified but `position` is omitted, then the default is
    `x` = b * BASE_WIDTH_SVG, `y` = [index of :any:`Helix`] * :any:`scadnano.distance_between_helices_svg`."""

    rotation: float = 0
    """Rotation angle (in degrees) of backbone of the :any:`Substrand` on this :any:`Helix` with 
    :py:data:`Substrand.forward` = ``True``. 
    
    The angle is relative to the offset :py:data:`Helix.rotation_anchor`, and 0 degrees is defined to
    be pointing *up* in both the side view and main view.
    
    A positive rotation angle rotates *clockwise* in the side view.
    This violates standard Cartesian coordinate conventions:
    https://en.wikipedia.org/wiki/Rotation_matrix, 
    but it is consistent with SVG rotation conventions:
    https://www.w3.org/TR/SVG11/coords.html#ExampleRotateScale.
    
    For example, a rotation of 90 degrees points right in the side view 
    and out of the screen in the main view.
    
    Default is 0 degrees."""

    rotation_anchor: int = 0
    """Offset on this :any:`Helix` that is the reference point for 0 degrees.
    The rotation at offset ``o`` is 360 degrees times the remainder of ``o - rotation_anchor`` 
    when divided by 10.5.
    
    For example, if :py:data:`Helix.rotation` = 0 and :py:data:`Helix.rotation_anchor` = 42, then
    at offsets of the form :math:`42 + 21k` for integer :math:`k` 
    (i.e., 42 itself, as well as 21, 0, -21, -42, ..., 63, 84, 105, ...),
    the rotation angle is also 0 at those offsets since
    they are integer multiples of 21 (hence also multiples of 10.5) from 42.
    
    Default is 0."""

    position3d: Position3D = None
    """Position (x,y,z) and orientation (pitch,roll,yaw) of this :any:`Helix` in 3D space.
    
    Optional if :py:data:`Helix.grid_position` is specified. 
    Default is pitch, roll, yaw are 0, and x,y,z are determined by grid position h, v, b."""

    _idx: int = -1

    # for optimization; list of substrands on that Helix
    _substrands: List['Substrand'] = field(default_factory=list)

    def to_json_serializable(self, suppress_indent=True):
        dct = dict()

        if self.min_offset != 0:
            dct[min_offset_key] = self.min_offset

        dct[max_offset_key] = self.max_offset

        if self.position3d is None:
            if self.grid_position[2] == 0:  # don't bother writing grid position base coordinate if it is 0
                dct[grid_position_key] = (self.grid_position[0], self.grid_position[1])
            else:
                dct[grid_position_key] = (self.grid_position[0], self.grid_position[1], self.grid_position[2])
        else:
            dct[position3d_key] = self.position3d.to_json_serializable(suppress_indent)

        # print(f'self.svg_position()    = {self.svg_position}')
        # print(f'default_svg_position() = {self.default_svg_position()}')
        default_x, default_y = self.default_svg_position()
        if not (_is_close(self.svg_position[0], default_x) and _is_close(self.svg_position[1], default_y)):
            dct[svg_position_key] = (self.svg_position[0], self.svg_position[1])

        if self.major_tick_distance is not None and self.major_tick_distance > 0:
            dct[major_tick_distance_key] = self.major_tick_distance

        if self.major_ticks is not None:
            dct[major_ticks_key] = self.major_ticks

        if self.rotation != 0:
            dct[rotation_key] = self.rotation

        if self.rotation_anchor != 0:
            dct[rotation_anchor_key] = self.rotation_anchor

        return _NoIndent(dct) if suppress_indent else dct

    def __post_init__(self):
        if self.major_ticks is not None:
            for major_tick in self.major_ticks:
                if major_tick > self.max_offset - self.min_offset:
                    raise IllegalDNADesignError(f'major tick {major_tick} in list {self.major_ticks} is '
                                                f'outside the range of available offsets since max_offset = '
                                                f'{self.max_offset}')

    def default_svg_position(self):
        return 0, self._idx * distance_between_helices_svg

    def set_idx(self, new_idx):
        self._idx = new_idx

    def idx(self):
        return self._idx

    def default_grid_position(self):
        return (0, self.idx(), 0)

    def calculate_major_ticks(self, default_major_tick_distance):
        """
        Calculates full list of major tick marks, whether using `default_major_tick_distance` (from
        :any:`DNADesign`), :py:data:`Helix.major_tick_distance`, or :py:data:`Helix.major_ticks`.
        They are used in reverse order to determine precedence. (e.g., :py:data:`Helix.major_ticks`
        overrides :py:data:`Helix.major_tick_distance`, which overrides
        `default_major_tick_distance` from :any:`DNADesign`.
        """
        if self.major_ticks is not None:
            return self.major_ticks
        distance = default_major_tick_distance if self.major_tick_distance <= 0 else self.major_tick_distance
        return list(range(self.min_offset, self.max_offset + 1, distance))

    @staticmethod
    def from_json(json_map: dict) -> Helix:
        grid_position = None
        if grid_position_key in json_map:
            gp_list = json_map[grid_position_key]
            if len(gp_list) not in [2, 3]:
                raise IllegalDNADesignError("list of grid_position coordinates must be length 2 or 3, "
                                            f"but this is the list: {gp_list}")
            if len(gp_list) == 2:
                gp_list.append(0)
            grid_position = tuple(gp_list)

        svg_position = None
        if svg_position_key in json_map:
            sp_list = json_map[grid_position_key]
            if len(sp_list) != 2:
                raise IllegalDNADesignError("svg_position must have exactly two integers, "
                                            f"but instead it has {len(sp_list)}: {sp_list}")
            svg_position = tuple(sp_list)

        major_tick_distance = json_map.get(major_tick_distance_key)
        major_ticks = json_map.get(major_ticks_key)
        min_offset = json_map.get(min_offset_key)
        max_offset = json_map.get(max_offset_key)
        rotation = json_map.get(rotation_key, default_helix_rotation)
        rotation_anchor = json_map.get(rotation_anchor_key, default_helix_rotation_anchor)
        position3d = json_map.get(position3d_key)

        return Helix(
            major_tick_distance=major_tick_distance,
            major_ticks=major_ticks,
            grid_position=grid_position,
            svg_position=svg_position,
            min_offset=min_offset,
            max_offset=max_offset,
            rotation=rotation,
            rotation_anchor=rotation_anchor,
            position3d=position3d,
        )


def _is_close(x1: float, x2: float):
    return abs(x1 - x2) < 0.00000001


@dataclass
class Substrand(_JSONSerializable):
    """
    A maximal portion of a :any:`Strand` that is continguous on a single :any:`Helix`.
    A :any:`Strand` contains a list of :any:`Substrand`'s (and also potentially :any:`Loopout`'s).
    """

    helix: int
    """index of the :any:`Helix` on which this :any:`Substrand` resides 
    in the list :any:`DNADesign.helices`."""

    forward: bool
    """Whether the strand "points" forward (i.e., its 3' end has a larger offset than its 5' end).
    If :any:`Substrand.forward` is ``True``, then 
    :any:`Substrand.start` is the 5' end of the :any:`Substrand` and 
    :any:`Substrand.end` is the 3' end of the :any:`Substrand`.
    If :any:`Substrand.forward` is ``False``, these roles are reversed."""

    start: int
    """
    The smallest offset position of any base on this Substrand
    (3' end if :any:`Substrand.forward` = ``False``,
    5' end if :any:`Substrand.forward` = ``True``).
    """

    # TODO: give option to user in constructor to specify that end is inclusive (default exclusive)
    end: int
    """
    1 plus the largest offset position of any base on this Substrand
    (5' end if :any:`Substrand.forward` = ``False``,
    3' end if :any:`Substrand.forward` = ``True``).
    Note that the set of base offsets occupied by this Substrand is {start, start+1, ..., end-1},
    i.e., inclusive for :py:data:`Strand.start` but exclusive for :py:data:`Strand.end`,
    the same convention used in Python for slices of lists and strings.
    (e.g., :samp:`"abcdef"[1:3] == "bc"`)
    
    Some methods (such as :py:meth:`Substrand.dna_sequence_in`) use the convention of being inclusive on 
    both ends and are marked with the word "INCLUSIVE".
    (Such a convention is easier to reason about when there are insertions and deletions.)
    """

    deletions: List[int] = field(default_factory=list)
    """List of positions of deletions on this Substrand."""

    insertions: List[Tuple[int, int]] = field(default_factory=list)
    """List of (position,num_insertions) pairs on this Substrand.
    
    This is the number of *extra* bases in addition to the base already at this position. 
    The total number of bases at this offset is num_insertions+1."""

    # not serialized; for efficiency
    _parent_strand: Strand = field(init=False, repr=False, compare=False, default=None)

    def __post_init__(self):
        self._check_start_end()

    def __repr__(self):
        rep = (f'Substrand(helix={self.helix}'
               f', forward={self.forward}'
               f', start={self.start}'
               f', end={self.end}') + \
              (f', deletions={self.deletions}' if len(self.deletions) > 0 else '') + \
              (f', insertions={self.insertions}' if len(self.insertions) > 0 else '') + \
              ')'
        return rep

    def __str__(self):
        return repr(self)

    def strand(self) -> Strand:
        return self._parent_strand

    def _check_start_end(self):
        if self.start >= self.end:
            raise StrandError(self._parent_strand,
                              f'start = {self.start} must be less than end = {self.end}')

    def to_json_serializable(self, suppress_indent=True):
        dct = OrderedDict()
        dct[helix_idx_key] = self.helix
        dct[forward_key] = self.forward
        dct[start_key] = self.start
        dct[end_key] = self.end
        if len(self.deletions) > 0:
            dct[deletions_key] = self.deletions
        if len(self.insertions) > 0:
            dct[insertions_key] = self.insertions
        return _NoIndent(dct) if suppress_indent else dct

    @staticmethod
    def is_loopout() -> bool:
        """Indicates if this is a :any:`Loopout` (always false)
        Useful when object could be either :any:`Loopout` or :any:`Substrand`."""
        return False

    @staticmethod
    def is_substrand() -> bool:
        """Indicates if this is a :any:`Substrand` (always true)
        Useful when object could be either :any:`Loopout` or :any:`Substrand`."""
        return True

    def set_start(self, new_start: int):
        self.start = new_start
        self._check_start_end()

    def set_end(self, new_end: int):
        self.end = new_end
        self._check_start_end()

    def offset_5p(self) -> int:
        """5' offset of this :any:`Substrand`, INCLUSIVE."""
        if self.forward:
            return self.start
        else:
            return self.end - 1
        # return self.start if self.forward else self.end - 1

    def offset_3p(self) -> int:
        """3' offset of this :any:`Substrand`, INCLUSIVE."""
        return self.end - 1 if self.forward else self.start

    def _num_insertions(self) -> int:
        # total number of insertions in this Substrand
        return sum(insertion[1] for insertion in self.insertions)

    def contains_offset(self, offset: int) -> bool:
        """Indicates if `offset` is the offset of a base on this substrand.

        Note that offsets refer to visual portions of the displayed grid for the Helix.
        If for example, this Substrand starts at position 0 and ends at 10, and it has 5 deletions,
        then it contains the offset 7 even though there is no base 7 positions from the start."""
        return self.start <= offset < self.end

    def __len__(self):
        """Same as :meth:`Substrand.dna_length`.

        See also :meth:`Substrand.visual_length`."""
        return self.dna_length()

    def dna_length(self) -> int:
        """Number of bases in this Substrand."""
        return self.end - self.start - len(self.deletions) + self._num_insertions()

    def dna_length_in(self, left, right) -> int:
        """Number of bases in this Substrand between `left` and `right` (INCLUSIVE)."""
        if not left <= right + 1:
            raise ValueError(f'left = {left} and right = {right} but we should have left <= right + 1')
        if not self.start <= left:
            raise ValueError(f'left = {left} should be at least self.start = {self.start}')
        if not right < self.end:
            raise ValueError(f'right = {right} should be at most self.end - 1 = {self.end - 1}')
        num_deletions = sum(1 for offset in self.deletions if left <= offset <= right)
        num_insertions = sum(length for (offset, length) in self.insertions if left <= offset <= right)
        return (right - left + 1) - num_deletions + num_insertions

    def visual_length(self) -> int:
        """Distance between :any:`Substrand.start` offset and :any:`Substrand.end` offset.

        This can be more or less than the :meth:`Substrand.dna_length` due to insertions and deletions."""
        return self.end - self.start

    def dna_sequence(self) -> Optional[str]:
        """Return DNA sequence of this Substrand, or ``None`` if no DNA sequence has been assigned
        to this :any:`Substrand`'s :any:`Strand`."""
        return self.dna_sequence_in(self.start, self.end - 1)

    def dna_sequence_in(self, offset_left: int, offset_right: int) -> Optional[str]:
        """Return DNA sequence of this Substrand in the interval of offsets given by
        [`offset_left`, `offset_right`], INCLUSIVE, or ``None`` if no DNA sequence has been assigned
        to this :any:`Substrand`'s :any:`Strand`.

        WARNING: This is inclusive on both ends,
        unlike other parts of this API where the right endpoint is exclusive.
        This is to make the notion well-defined when one of the endpoints is on an offset with a
        deletion or insertion."""
        strand_seq = self._parent_strand.dna_sequence
        if strand_seq is None:
            return None

        # if on a deletion, move inward until we are off of it
        while offset_left in self.deletions:
            offset_left += 1
        while offset_right in self.deletions:
            offset_right -= 1

        if offset_left > offset_right:
            return ''
        if offset_left >= self.end:
            return ''
        if offset_right < 0:
            return ''

        str_idx_left = self.offset_to_str_idx(offset_left, self.forward)
        str_idx_right = self.offset_to_str_idx(offset_right, not self.forward)
        if not self.forward:  # these will be out of order if strand is left
            str_idx_left, str_idx_right = str_idx_right, str_idx_left
        subseq = strand_seq[str_idx_left:str_idx_right + 1]
        return subseq

    def get_seq_start_idx(self) -> int:
        """Starting DNA subsequence index for first base of this :any:`Substrand` on its
        Parent :any:`Strand`'s DNA sequence."""
        substrands = self._parent_strand.substrands
        # index of self in parent strand's list of substrands
        self_substrand_idx = substrands.index(self)
        # index of self's position within the DNA sequence of parent strand
        self_seq_idx_start = sum(prev_substrand.dna_length()
                                 for prev_substrand in substrands[:self_substrand_idx])
        return self_seq_idx_start

    def offset_to_str_idx(self, offset: int, offset_closer_to_5p: bool) -> int:
        """ Convert from offset on this :any:`Substrand`'s :any:`Helix`
        to string index on the parent :any:`Strand`'s DNA sequence.

        If `offset_closer_to_5p` is ``True``, (this only matters if `offset` contains an insertion)
        then the only leftmost string index corresponding to this offset is included,
        otherwise up to the rightmost string index (including all insertions) is included."""
        if offset in self.deletions:
            raise ValueError(f'offset {offset} illegally contains a deletion from {self.deletions}')

        # length adjustment for insertions depends on whether this is a left or right offset
        len_adjust = self._net_ins_del_length_increase_from_5p_to(offset, offset_closer_to_5p)

        # get string index assuming this Substrand is first on Strand
        if self.forward:
            offset += len_adjust  # account for insertions and deletions
            ss_str_idx = offset - self.start
        else:
            # account for insertions and deletions
            offset -= len_adjust  # account for insertions and deletions
            ss_str_idx = self.end - 1 - offset

        # correct for existence of previous Substrands on this Strand
        return ss_str_idx + self.get_seq_start_idx()

    def _net_ins_del_length_increase_from_5p_to(self, offset_edge: int, offset_closer_to_5p: bool) -> int:
        """Net number of insertions from 5'/3' end to offset_edge,
        INCLUSIVE on 5'/3' end, EXCLUSIVE on offset_edge.

        Set `five_p` ``= False`` to test from 3' end to `offset_edge`."""
        length_increase = 0
        for deletion in self.deletions:
            if self._between_5p_and_offset(deletion, offset_edge):
                length_increase -= 1
        for (insertion_offset, insertion_length) in self.insertions:
            if self._between_5p_and_offset(insertion_offset, offset_edge):
                length_increase += insertion_length
        # special case for when offset_edge is an endpoint closer to the 3' end,
        # we add its extra insertions also in this case
        if not offset_closer_to_5p:
            insertion_map: Dict[int, int] = dict(self.insertions)
            if offset_edge in insertion_map:
                insertion_length = insertion_map[offset_edge]
                length_increase += insertion_length
        return length_increase

    def _between_5p_and_offset(self, offset_to_test: int, offset_edge: int) -> bool:
        return ((self.forward and self.start <= offset_to_test < offset_edge) or
                (not self.forward and offset_edge < offset_to_test < self.end))

    # def _between_3p_and_offset(self, offset_to_test: int, offset_edge: int) -> bool:
    #     return ((self.direction == Direction.left and self.start <= offset_to_test < offset_edge) or
    #             (self.direction == Direction.forward and offset_edge < offset_to_test < self.end))

    # The type hint 'Substrand' must be in quotes since Substrand is not yet defined.
    # This is a "forward reference": https://www.python.org/dev/peps/pep-0484/#forward-references
    def overlaps(self, other: Substrand) -> bool:
        r"""Indicates if this substrand's set of offsets (the set
        :math:`\{x \in \mathbb{N} \mid`
        ``self.start``
        :math:`\leq x \leq`
        ``self.end``
        :math:`\}`)
        has nonempty intersection with those of `other`,
        and they appear on the same helix,
        and they point in opposite directions."""  # noqa (suppress PEP warning)
        return (self.helix == other.helix and
                self.forward == (not other.forward) and
                self.compute_overlap(other)[0] >= 0)

    def overlaps_illegally(self, other: Substrand):
        r"""Indicates if this substrand's set of offsets (the set
        :math:`\{x \in \mathbb{N} \mid`
        ``self.start``
        :math:`\leq x \leq`
        ``self.end``
        :math:`\}`)
        has nonempty intersection with those of `other`,
        and they appear on the same helix,
        and they point in the same direction."""  # noqa (suppress PEP warning)
        return (self.helix == other.helix and
                self.forward == other.forward and
                self.compute_overlap(other)[0] >= 0)

    def compute_overlap(self, other: Substrand) -> Tuple[int, int]:
        """Return [left,right) offset indicating overlap between this Substrand and `other`.

        Return ``(-1,-1)`` if they do not overlap (different helices, or non-overlapping regions
        of the same helix)."""
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        if overlap_start >= overlap_end:  # overlap is empty
            return -1, -1
        return overlap_start, overlap_end

    def insertion_offsets(self) -> List[int]:
        """Return offsets of insertions (but not their lengths)."""
        return [ins_off for (ins_off, _) in self.insertions]

    @staticmethod
    def from_json(json_map):
        helix = json_map[helix_idx_key]
        forward = json_map[forward_key]
        start = json_map[start_key]
        end = json_map[end_key]
        deletions = json_map.get(deletions_key, [])
        insertions = list(map(tuple, json_map.get(insertions_key, [])))
        return Substrand(
            helix=helix,
            forward=forward,
            start=start,
            end=end,
            deletions=deletions,
            insertions=insertions,
        )


'''
    var forward = util.get_value(json_map, constants.forward_key, name);
    var helix = util.get_value(json_map, constants.helix_idx_key, name);
    var start = util.get_value(json_map, constants.start_key, name);
    var end = util.get_value(json_map, constants.end_key, name);
//    List<int> deletions =
//        json_map.containsKey(constants.deletions_key) ? List<int>.from(json_map[constants.deletions_key]) : [];
//    List<Tuple2<int, int>> insertions =
//        json_map.containsKey(constants.insertions_key) ? parse_json_insertions(json_map[constants.insertions_key]) : [];
    var deletions = List<int>.from(util.get_value_with_default(json_map, constants.deletions_key, []));
    var insertions =
        parse_json_insertions(util.get_value_with_default(json_map, constants.insertions_key, []));
'''


@dataclass
class Loopout(_JSONSerializable):
    """Represents a single-stranded loopout on a :any:`Strand`.

    One could think of a :any:`Loopout` as a type of :any:`Substrand`, but none of the fields of
    :any:`Substrand` make sense for :any:`Loopout`, so they are not related to each other in the type
    hierarchy. It is interpreted that a :any:`Loopout` is a single-stranded region bridging two
    :any:`Substrand`'s that are connected to :any:`Helix`'s, or if it occurs on the end of a :any:`Strand`,
    then it is a single-stranded extension. It is illegal for two consecutive :any:`Substrand`'s to both
    be :any:`Loopout`'s, and for a :any:`Strand` to have only one element of :any:`Strand.substrands`
    that is a :any:`Loopout`.

    Loopout has only a single field :py:data:`Loopout.length` that specifies the length of the loopout.

    For example, one use of a loopout is to describe a hairpin (a.k.a.,
    `stem-loop <https://en.wikipedia.org/wiki/Stem-loop>`_).
    The following creates a :any:`Strand` that represents a hairpin with a stem length of 10 and a loop
    length of 5.

    .. code-block:: Python

        import scadnano as sc

        ss_f = sc.Substrand(helix=0, forward=True, start=0, end=10)
        loop = sc.Loopout(length=5)
        ss_r = sc.Substrand(helix=0, forward=False, start=0, end=10)
        hairpin = sc.Strand([ss_f, loop, ss_r])
    """

    length: int
    """Length (in DNA bases) of this Loopout."""

    # not serialized; for efficiency
    _parent_strand: Strand = field(init=False, repr=False, compare=False, default=None)

    def to_json_serializable(self, suppress_indent=True):
        dct = {loopout_key: self.length}
        return _NoIndent(dct)

    def __repr__(self):
        return f'Loopout({self.length})'

    def __str__(self):
        return repr(self)

    @staticmethod
    def is_loopout() -> bool:
        """Indicates if this is a :any:`Loopout` (always true).
        Useful when object could be either :any:`Loopout` or :any:`Substrand`."""
        return True

    @staticmethod
    def is_substrand() -> bool:
        """Indicates if this is a :any:`Substrand` (always false)
        Useful when object could be either :any:`Loopout` or :any:`Substrand`."""
        return False

    def __len__(self):
        """Same as :any:`Loopout.dna_length`"""
        return self.dna_length()

    def dna_length(self) -> int:
        """Length of this :any:`Loopout`; same as field :py:data:`Loopout.length`."""
        return self.length

    def dna_sequence(self) -> Optional[str]:
        """Return DNA sequence of this :any:`Loopout`, or ``None`` if no DNA sequence has been assigned
        to the :any:`Strand` of this :any:`Loopout`."""
        strand_seq = self._parent_strand.dna_sequence
        if strand_seq is None:
            return None

        str_idx_left = self.get_seq_start_idx()
        str_idx_right = str_idx_left + self.length  # EXCLUSIVE (unlike similar code for Substrand)
        subseq = strand_seq[str_idx_left:str_idx_right]
        return subseq

    def get_seq_start_idx(self) -> int:
        """Starting DNA subsequence index for first base of this :any:`Loopout` on its
        :any:`Strand`'s DNA sequence."""
        substrands = self._parent_strand.substrands
        # index of self in parent strand's list of substrands
        self_substrand_idx = substrands.index(self)
        # index of self's position within the DNA sequence of parent strand
        self_seq_idx_start = sum(prev_substrand.dna_length()
                                 for prev_substrand in substrands[:self_substrand_idx])
        return self_seq_idx_start

    @staticmethod
    def from_json(json_map):
        if loopout_key not in json_map:
            raise IllegalDNADesignError(f'no key "{loopout_key}" in JSON map')
        length = int(json_map[loopout_key])
        return Loopout(length=length)


_wctable = str.maketrans('ACGTacgt', 'TGCAtgca')


def wc(seq: str) -> str:
    """Return reverse Watson-Crick complement of `seq`."""
    return seq.translate(_wctable)[::-1]


@dataclass
class IDTFields(_JSONSerializable):
    """Data required when ordering DNA strands from the synthesis company
    `IDT DNA Technologies <https://www.idtdna.com/>`_.
    This data is used when automatically generating files used to order DNA from IDT."""

    name: str
    """Name of the strand (first field in IDT bulk input: https://www.idtdna.com/site/order/oligoentry).
    
    Non-optional field.
    """

    scale: str = default_idt_scale
    """Synthesis scale at which to synthesize the strand (third field in IDT bulk input:
    https://www.idtdna.com/site/order/oligoentry).
    Choices supplied by IDT at the time this was written: 
    ``"25nm"``, ``"100nm"``, ``"250nm"``, ``"1um"``, ``"5um"``, 
    ``"10um"``, ``"4nmU"``, ``"20nmU"``, ``"PU"``, ``"25nmS"``.
    
    Optional field.
    """

    purification: str = default_idt_purification
    """Purification options (fourth field in IDT bulk input:
    https://www.idtdna.com/site/order/oligoentry). 
    Choices supplied by IDT at the time this was written: 
    ``"STD"``, ``"PAGE"``, ``"HPLC"``, ``"IEHPLC"``, ``"RNASE"``, ``"DUALHPLC"``, ``"PAGEHPLC"``.
    
    Optional field.
    """

    plate: Optional[str] = None
    """Name of plate in case this strand will be ordered on a 96-well or 384-well plate.
    
    Optional field, but non-optional if :py:data:`IDTField.well` is not ``None``.
    """

    well: Optional[str] = None
    """Well position on plate in case this strand will be ordered on a 96-well or 384-well plate.
    
    Optional field, but non-optional if :py:data:`IDTField.plate` is not ``None``.
    """

    def __post_init__(self):
        _check_idt_string_not_none_or_empty(self.name, 'name')
        _check_idt_string_not_none_or_empty(self.scale, 'scale')
        _check_idt_string_not_none_or_empty(self.purification, 'purification')
        if self.plate is None and self.well is not None:
            raise IllegalDNADesignError(f'IDTFields.plate cannot be None if IDTFields.well is not None\n'
                                        f'IDTFields.well = {self.well}')
        if self.plate is not None and self.well is None:
            raise IllegalDNADesignError(f'IDTFields.well cannot be None if IDTFields.plate is not None\n'
                                        f'IDTFields.plate = {self.plate}')

    def to_json_serializable(self, suppress_indent=True):
        dct = self.__dict__
        if self.plate is None:
            del dct['plate']
        if self.well is None:
            del dct['well']
        return _NoIndent(dct)


def _check_idt_string_not_none_or_empty(value: str, field_name: str):
    if value is None:
        raise IllegalDNADesignError(f'field {field_name} in IDTFields cannot be None')
    if len(value) == 0:
        raise IllegalDNADesignError(f'field {field_name} in IDTFields cannot be empty')


@dataclass
class Strand(_JSONSerializable):
    """
    Represents a single strand of DNA.

    Each maximal portion that is continguous on a single :any:`Helix` is a :any:`Substrand`.
    Crossovers from one :any:`Helix` to another are implicitly from the 3' end of one of this
    Strand's :any:`Substrand`'s to the 5' end of the next :any:`Substrand`.

    A portion of the :any:`Strand` not associated to any :any:`Helix` is represented by a :any:`Loopout`.
    Two :any:`Loopout`'s cannot occur consecutively on a :any:`Strand`, nor can a :any:`Strand`
    contain only a :any:`Loopout` but no :any:`Substrand`.

    To give a strand the same color that
    `cadnano <https://cadnano.org/>`_
    uses for the scaffold,
    use :any:`scadnano.default_scaffold_color` in the :any:`Strand` constructor:

    .. code-block:: Python

        import scadnano as sc

        scaffold_substrands = [ ... ]
        scaffold_strand = sc.Strand(substrands=scaffold_substrands, color=sc.default_scaffold_color)

    Alternately, one can use :any:`DNAOrigamiDesign`, which is a subclass of
    :any:`DNADesign` that has one extra field representing the scaffold strand.
    This class will automatically assign the color of the scaffold strand to be
    :py:data:`default_scaffold_color`.
    """

    substrands: List[Union[Substrand, Loopout]]
    """:any:`Substrand`'s (or :any:`Loopout`'s) composing this Strand. 
    Each :any:`Substrand` is contiguous on a single :any:`Helix` 
    and could be either single-stranded or double-stranded, 
    whereas each :any:`Loopout` is single-stranded and has no associated :any:`Helix`."""

    dna_sequence: Optional[str] = None
    """Do not assign directly to this field. Always use :any:`DNADesign.assign_dna` 
    (for complementarity checking) or :any:`Strand.set_dna_sequence` 
    (without complementarity checking, to allow mismatches)."""

    color: Optional[Color] = None
    """Color to show this strand in the main view. If not specified in the constructor,
    a color is assigned by cycling through a list of defaults given by 
    :meth:`ColorCycler.colors`"""

    automatically_assign_color: bool = field(repr=False, default=True)
    """If `automatically_assign_color` = ``False`` and `color` = ``None``, do not automatically
    assign a :any:`Color` to this :any:`Strand`. 
    In this case color will be set to its default of ``None`` and will not be
    written to the JSON with :py:meth:`DNADesign.write_scadnano_file` or :py:meth:`DNADesign.to_json`."""

    idt: Optional[IDTFields] = None
    """Fields used when ordering strands from the synthesis company IDT 
    (Integrated DNA Technologies, Coralville, IA). If present (i.e., not equal to :const:`None`)
    then the method :py:meth:`DNADesign.write_idt_bulk_input_file` can be called to automatically
    generate an text file for ordering strands in test tubes: 
    https://eu.idtdna.com/site/order/oligoentry,
    as can the method :py:meth:`DNADesign.write_idt_plate_excel_file` for writing a Microsoft Excel 
    file that can be uploaded to IDT's website for describing DNA sequences to be ordered in 96-well
    or 384-well plates."""

    use_default_idt: bool = False
    """If ``True``, assigns an :any:`IDTFields` to this :any:`Strand` with same naming convention as
    cadnano, i.e., :py:data:`IDTFields.name` = "ST{h5}[{s}]{h3}[{e}]", where h5 and h3 are the 
    :any:`Helix`'s of the 5' and 3' ends, respectively, of the :any:`Strand`, 
    and s and e are the respective start and end offsets on those helices.
    """

    # not serialized; efficient way to see a list of all substrands on a given helix
    _helix_idx_substrand_map: Dict[int, List[Substrand]] = field(
        init=False, repr=False, compare=False, default=None)

    def to_json_serializable(self, suppress_indent=True):
        dct = OrderedDict()
        if self.color is not None:
            dct[color_key] = self.color.to_json_serializable(suppress_indent)
        if self.dna_sequence is not None:
            dct[dna_sequence_key] = self.dna_sequence
        if self.idt is not None:
            dct[idt_key] = self.idt.to_json_serializable(suppress_indent)
        dct[substrands_key] = [substrand.to_json_serializable(suppress_indent) for substrand in
                               self.substrands]
        if hasattr(self, is_scaffold_key):
            dct[is_scaffold_key] = self.is_scaffold
        return dct

    def __post_init__(self):
        # if color not specified, pick one by cycling through list of staple colors,
        # unless caller specified not to
        global color_cycler
        if self.color is None and self.automatically_assign_color:
            self.color = next(color_cycler)

        self._helix_idx_substrand_map = defaultdict(list)

        for substrand in self.substrands:
            if substrand.is_substrand():
                self._helix_idx_substrand_map[substrand.helix].append(substrand)

        for substrand in self.substrands:
            substrand._parent_strand = self

        if len(self.substrands) == 1:
            if self.first_substrand().is_loopout():
                raise StrandError(self, 'strand cannot have a single Loopout as its only substrand')

        for ss1, ss2 in _pairwise(self.substrands):
            if ss1.is_loopout() and ss2.is_loopout():
                raise StrandError(self, 'cannot have two consecutive Loopouts in a strand')

        if self.use_default_idt:
            self.set_default_idt(True)

    def __eq__(self, other: Strand) -> bool:
        if not isinstance(other, Strand):
            return False
        return self.substrands == other.substrands

    def __hash__(self):
        return hash(self.substrands)

    def set_color(self, color: Color):
        """Sets color of this :any:`Strand`."""
        self.color = color;

    def set_default_idt(self, use_default_idt):
        """Sets idt field to be the default given the Substrand data of this :any:`Strand`."""
        self.use_default_idt = use_default_idt
        if use_default_idt:
            start_helix = self.first_bound_substrand().helix
            end_helix = self.last_bound_substrand().helix
            start_offset = self.first_bound_substrand().offset_5p()
            end_offset = self.last_bound_substrand().offset_3p()
            self.idt = IDTFields(name=f'ST{start_helix}[{start_offset}]{end_helix}[{end_offset}]')
        else:
            self.idt = None

    def first_substrand(self) -> Union[Substrand, Loopout]:
        """First substrand (of type either :any:`Substrand` or :any:`Loopout`) on this :any:`Strand`."""
        return self.substrands[0]

    def last_substrand(self) -> Union[Substrand, Loopout]:
        """Last substrand (of type either :any:`Substrand` or :any:`Loopout`) on this :any:`Strand`."""
        return self.substrands[-1]

    def set_dna_sequence(self, sequence: str):
        """Set this :any:`Strand`'s DNA sequence to `seq`
        WITHOUT checking for complementarity with overlapping
        :any:`Strand`'s or automatically assigning their sequences.
        To assign a sequence to a :any:`Strand` and have the overlapping
        :any:`Strand`'s automatically have the appropriate Watson-Crick complements assigned,
        use :any:`DNADesign.assign_dna`.

        All whitespace in `sequence` is removed,
        and lowercase bases 'a', 'c', 'g', 't' are converted to uppercase.

        `sequence`, after all whitespace is removed, must be exactly the same length as
        :py:meth:`Strand.dna_length`.
        Wildcard symbols (:py:const:`DNA_case_wildcard`) are allowed to leave part of the DNA unassigned.
        """
        trimmed_seq = _remove_whitespace_and_uppercase(sequence)
        if len(trimmed_seq) != self.dna_length():
            ss = self.first_substrand()
            raise StrandError(self, f"strand starting at helix {ss.helix} offset {ss.offset_5p()} "
                                    f"has length {self.dna_length()}, but you attempted to assign a "
                                    f"DNA sequence of length {len(trimmed_seq)}: {sequence}")
        self.dna_sequence = trimmed_seq

    def dna_length(self) -> int:
        """Return sum of DNA length of :any:`Substrand`'s and :any:`Loopout`'s of this :any:`Strand`."""
        acc = 0
        for substrand in self.substrands:
            acc += substrand.dna_length()
        return acc

    def bound_substrands(self) -> List[Substrand]:
        """:any:`Substrand`'s of this :any:`Strand` that are not :any:`Loopout`'s."""
        return [ss for ss in self.substrands if ss.is_substrand()]

    def offset_5p(self) -> int:
        """5' offset of this entire :any:`Strand`, INCLUSIVE."""
        return self.first_substrand().offset_5p()

    def offset_3p(self) -> int:
        """3' offset of this entire :any:`Strand`, INCLUSIVE."""
        return self.last_substrand().offset_3p()

    def overlaps(self, other: Strand) -> bool:
        """Indicates whether `self` overlaps `other_strand`, meaning that the set of offsets occupied
        by `self` has nonempty intersection with those occupied by `other_strand`."""
        for substrand_self in self.bound_substrands():
            for substrand_other in other.bound_substrands():
                if substrand_self.overlaps(substrand_other):
                    return True
        return False

    def assign_dna_complement_from(self, other: Strand):
        """Assuming a DNA sequence has been assigned to `other`, assign its Watson-Crick
        complement to the portions of this Strand that are bound to `other`.

        Generally this is not called directly; use :py:meth:`DNADesign.assign_dna` to assign
        a DNA sequence to a :any:`Strand`. The method :py:meth:`DNADesign.assign_dna` will calculate
        which other :any:`Strand`'s need
        to be assigned via :py:meth:`Strand.assign_dna_complement_from`.

        However, it is permitted to assign the field :py:data:`Strand.dna_sequence` directly
        via the method :py:meth:`Strand.set_dna_sequence`.
        This is used, for instance, to assign a DNA sequence to a :any:`Strand` bound to another
        :any:`Strand`
        with an assigned DNA sequence where they overlap. In this case no error checking
        about sequence complementarity is done. This can be used to intentionally assign *mismatching*
        DNA sequences to :any:`Strand`'s that are bound on a :any:`Helix`."""

        already_assigned = self.dna_sequence is not None

        # put DNA sequences to assign to substrands in List, one position per substrand
        strand_complement_builder = []
        if already_assigned:
            for substrand in self.substrands:
                strand_complement_builder.append(substrand.dna_sequence())
        else:
            for substrand in self.substrands:
                wildcards = DNA_base_wildcard * substrand.dna_length()
                strand_complement_builder.append(wildcards)

        for (ss_idx, substrand_self) in enumerate(self.substrands):
            if substrand_self.is_loopout():
                substrand_self_dna_sequence = DNA_base_wildcard * substrand_self.dna_length()
            else:
                helix = substrand_self.helix

                # for helix, substrands_on_helix_self in self._helix_idx_substrand_map.items():
                substrands_on_helix_other = other._helix_idx_substrand_map[helix]
                # for substrand_self in substrands_on_helix_self:
                overlaps = []
                for substrand_other in substrands_on_helix_other:
                    if substrand_self != substrand_other and substrand_self.overlaps(substrand_other):
                        overlap = substrand_self.compute_overlap(substrand_other)
                        overlaps.append((overlap, substrand_other))

                overlaps.sort()

                substrand_complement_builder = []
                start_idx = substrand_self.start
                # repeatedly insert wildcards into gaps, then reverse WC complement
                for ((overlap_left, overlap_right), substrand_other) in overlaps:
                    # wildcards = DNA_base_wildcard * (overlap_left - start_idx)
                    num_wildcard_bases = substrand_self.dna_length_in(start_idx, overlap_left - 1)
                    wildcards = DNA_base_wildcard * num_wildcard_bases

                    other_seq = substrand_other.dna_sequence_in(overlap_left, overlap_right - 1)
                    overlap_complement = wc(other_seq)
                    substrand_complement_builder.append(wildcards)
                    substrand_complement_builder.append(overlap_complement)
                    start_idx = overlap_right

                # last wildcard for gap between last overlap and end
                # last_wildcards = DNA_base_wildcard * (substrand_self.end - start_idx)
                num_wildcard_bases = substrand_self.dna_length_in(start_idx, substrand_self.end - 1)
                last_wildcards = DNA_base_wildcard * num_wildcard_bases

                substrand_complement_builder.append(last_wildcards)

                # If pointing left, each individual overlap sequence was reverse orientation in wc(),
                # but not the list of all of them put together until now.
                if not substrand_self.forward:
                    substrand_complement_builder.reverse()

                substrand_self_dna_sequence = ''.join(substrand_complement_builder)

            # merge with existing pre-assigned sequence
            existing_substrand_self_dna_sequence = strand_complement_builder[ss_idx]
            merged_substrand_self_dna_sequence = _string_merge_wildcard(substrand_self_dna_sequence,
                                                                        existing_substrand_self_dna_sequence,
                                                                        DNA_base_wildcard)
            strand_complement_builder[ss_idx] = merged_substrand_self_dna_sequence

        strand_complement = ''.join(strand_complement_builder)
        new_dna_sequence = strand_complement
        if self.dna_sequence is not None:
            try:
                new_dna_sequence = _string_merge_wildcard(self.dna_sequence, new_dna_sequence,
                                                          DNA_base_wildcard)
            except ValueError:
                ss_self = self.first_substrand()
                ss_other = other.first_substrand()
                msg = f'strand starting at helix {ss_self.helix}, offset {ss_self.offset_5p()} has ' \
                      f'length ' \
                      f'{self.dna_length()} and already has a partial DNA sequence assignment of length ' \
                      f'{len(self.dna_sequence)}, which is \n' \
                      f'{self.dna_sequence}, ' \
                      f'but you tried to assign sequence of length {len(new_dna_sequence)} to it, which ' \
                      f'is\n{new_dna_sequence} (this assignment was indirect, since you assigned directly ' \
                      f'to a strand bound to this one). This occurred while directly assigning a DNA ' \
                      f'sequence to the strand whose 5\' end is at helix {ss_other.helix}, and is of ' \
                      f'length {other.dna_length()}.'
                raise IllegalDNADesignError(msg)

        self.set_dna_sequence(new_dna_sequence)
        # self.dna_sequence = _pad_dna(new_dna_sequence, self.dna_length())

    def _insert_substrand(self, order, substrand):
        """Only intended to be called by DNADesign.insert_substrand"""
        self.substrands.insert(order, substrand)
        substrand._parent_strand = self
        if substrand.is_substrand():
            self._helix_idx_substrand_map[substrand.helix].append(substrand)
        if self.use_default_idt:
            self.set_default_idt()

    def _remove_substrand(self, substrand):
        """Only intended to be called by DNADesign.remove_substrand"""
        self.substrands.remove(substrand)
        substrand._parent_strand = None
        if substrand.is_substrand():
            self._helix_idx_substrand_map[substrand.helix].remove(substrand)
        if self.use_default_idt:
            self.set_default_idt()

    def contains_loopouts(self):
        for ss in self.substrands:
            if ss.is_loopout():
                return True
        return False

    def first_bound_substrand(self):
        """First :any:`Substrand` (i.e., not a :any:`Loopout`) on this :any:`Strand`.

        Currently the first and last strand must not be :any:`Loopout`'s, so this should return the same
        substrand as :py:meth:`Strand.first_substrand`, but in case an initial or final :any:`Loopout` is
        supported in the future, this method is provided."""
        for substrand in self.substrands:
            if substrand.is_substrand():
                return substrand

    def last_bound_substrand(self):
        """Last :any:`Substrand` (i.e., not a :any:`Loopout`) on this :any:`Strand`.

        Currently the first and last strand must not be :any:`Loopout`'s, so this should return the same
        substrand as :py:meth:`Strand.first_substrand`, but in case an initial or final :any:`Loopout` is
        supported in the future, this method is provided."""
        substrands_rev = list(self.substrands)
        substrands_rev.reverse()
        for substrand in substrands_rev:
            if substrand.is_substrand():
                return substrand

    def reverse(self):
        """
        Reverses "polarity" of this :any:`Strand`.

        Does NOT check whether this keeps the :any:`DNADesign` legal, so be cautious in calling this method
        directly. To reverse every :any:`Strand`, called :py:meth:`DNADesign.reverse_all`.
        If the design was legal before, it will be legal after calling that method.
        """
        self.substrands.reverse()
        for substrand in self.bound_substrands():
            substrand.forward = not substrand.forward

    @staticmethod
    def from_json(json_map: dict) -> Strand:
        if substrands_key not in json_map:
            raise IllegalDNADesignError(f'key "{substrands_key}" is missing from the description of a Strand:'
                                        f'\n  {json_map}')
        substrand_jsons = json_map[substrands_key]
        if len(substrand_jsons) == 0:
            raise IllegalDNADesignError(f'substrands list cannot be empty')

        substrands = []
        for substrand_json in substrand_jsons:
            if loopout_key in substrand_json:
                substrands.append(Loopout.from_json(substrand_json))
            else:
                substrands.append(Substrand.from_json(substrand_json))
        if isinstance(substrands[0], Loopout):
            raise IllegalDNADesignError('Loopout at beginning of Strand not supported')
        if isinstance(substrands[-1], Loopout):
            raise IllegalDNADesignError('Loopout at end of Strand not supported')

        is_scaffold = json_map.get(is_scaffold_key, False)
        dna_sequence = json_map.get(dna_sequence_key)
        idt = json_map.get(idt_key)
        color_str = json_map.get(color_key, default_scaffold_color if is_scaffold else default_strand_color)
        color = Color(hex=color_str)

        return Strand(
            substrands=substrands,
            dna_sequence=dna_sequence,
            color=color,
            idt=idt,
        )


def _string_merge_wildcard(s1: str, s2: str, wildcard: str) -> str:
    """Takes a "union" of two equal-length strings `s1` and `s2`.
    Whenever one has a symbol `wildcard` and the other does not, the result has the non-wildcard symbol.

    Raises :py:class:`ValueError` if `s1` and `s2` are not the same length or do not agree on non-wildcard
    symbols at any position."""
    if len(s1) != len(s2):
        raise ValueError(f'\ns1={s1} and\ns2={s2}\nare not the same length.')
    union_builder = []
    for i in range(len(s1)):
        c1, c2 = s1[i], s2[i]
        if c1 == wildcard:
            union_builder.append(c2)
        elif c2 == wildcard:
            union_builder.append(c1)
        elif c1 != c2:
            raise ValueError(f's1={s1} and s2={s2} have unequal symbols {c1} and {c2} at position {i}.')
        elif c1 == c2:
            union_builder.append(c1)
        else:
            raise AssertionError('should be unreachable')
    return ''.join(union_builder)


class IllegalDNADesignError(ValueError):
    """Indicates that some aspect of the :any:`DNADesign` object is illegal."""

    def __init__(self, the_cause: str):
        self.cause = the_cause

    # __str__ is to print() the value
    def __str__(self):
        return repr(self.cause)


class StrandError(IllegalDNADesignError):
    """Indicates that the :any:`DNADesign` is illegal due to some specific :any:`Strand`.
    Information about the :any:`Strand` is embedded in the error message when this exception is
    raised that helps to identify which :any:`Strand` caused the problem."""

    def __init__(self, strand: Strand, the_cause: str):
        first_substrand = strand.first_bound_substrand()
        last_substrand = strand.last_bound_substrand()

        msg = (f'{the_cause}\n'
               f'strand length        =  {strand.dna_length()}\n'
               f'DNA length           =  {len(strand.dna_sequence) if strand.dna_sequence else "N/A"}\n'
               f'DNA sequence         =  {strand.dna_sequence}'
               f"strand 5' helix      =  {first_substrand.helix if first_substrand else 'N/A'}\n"
               f"strand 5' end offset =  {first_substrand.offset_5p() if first_substrand else 'N/A'}\n"
               f"strand 3' helix      =  {last_substrand.helix if last_substrand else 'N/A'}\n"
               f"strand 3' end offset =  {last_substrand.offset_3p() if last_substrand else 'N/A'}\n")

        super().__init__(msg)
        # super(IllegalDNADesignError, self).__init__(msg)


# TODO: add mutation operations to DNADesign to mutate all of its parts:
#  - Helix
#    - idx
#    - max_bases
#    - major_ticks (possibly as part of a deletion on strands that actually deletes the whole offset)
#    - svg_position (can help with importing cadnano designs that display poorly with default positions)
#  - Substrand
#    - helix
#    - right
#    - start
#    - end
#  - Strand
#    - unassign DNA sequence
#  - DNADesign
#    - add Helix
#    - remove Helix


def _plates(idt_strands):
    plates = set()
    for strand in idt_strands:
        if strand.idt is not None and strand.idt.plate is not None:
            plates.add(strand.idt.plate)
    return list(plates)


_96WELL_PLATE_ROWS: List[str] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
_96WELL_PLATE_COLS: List[int] = list(range(1, 13))

_384WELL_PLATE_ROWS: List[str] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                  'P']
_384WELL_PLATE_COLS: List[int] = list(range(1, 25))


@enum.unique
class PlateType(int, enum.Enum):
    """Represents two different types of plates in which DNA sequences can be ordered."""

    wells96 = 96
    """96-well plate."""

    wells384 = 384
    """384-well plate."""

    def rows(self) -> List[str]:
        return _96WELL_PLATE_ROWS if self is PlateType.wells96 else _384WELL_PLATE_ROWS

    def cols(self) -> List[int]:
        return _96WELL_PLATE_COLS if self is PlateType.wells96 else _384WELL_PLATE_COLS


class _PlateCoordinate:

    def __init__(self, plate_type: PlateType):
        self._plate_type = plate_type
        self._plate: int = 1
        self._row_idx: int = 0
        self._col_idx: int = 0

    def increment(self):
        self._row_idx += 1
        if self._row_idx == len(self._plate_type.rows()):
            self._row_idx = 0
            self._col_idx += 1
            if self._col_idx == len(self._plate_type.cols()):
                self._col_idx = 0
                self._plate += 1

    def plate(self) -> int:
        return self._plate

    def row(self) -> str:
        return self._plate_type.rows()[self._row_idx]

    def col(self) -> int:
        return self._plate_type.cols()[self._col_idx]

    def well(self) -> str:
        return f'{self.row()}{self.col()}'


@dataclass
class DNADesign(_JSONSerializable):
    """Object representing the entire design of the DNA structure."""

    strands: List[Strand]
    """All of the :any:`Strand`'s in this :any:`DNADesign`.
    
    Required field."""

    helices: List[Helix] = None
    """All of the :any:`Helix`'s in this :any:`DNADesign`. 
    
    Optional field. If not specified, then the number of helices will be just large enough to store the
    largest index :py:data:`Substrand.helix` 
    stored in any :any:`Substrand` 
    in :py:data:`DNADesign.strands`."""

    grid: Grid = Grid.square
    """Common choices for how to arrange helices relative to each other.
    
    Optional field."""

    major_tick_distance: int = -1
    """Distance between major ticks (bold) delimiting boundaries between bases.
    
    Optional field.
    If not specified, default value is 8 unless overridden by :py:data:`DNADesign.grid`.
    If 0 then no major ticks are drawn.
    If negative then the default value is assumed, but `major_tick_distance` is not stored in the JSON file
    when serialized.
    If :any:`DNADesign.grid` = :any:`Grid.square` then the default value is 8.
    If :any:`DNADesign.grid` = :any:`Grid.hex` or :any:`Grid.honeycomb` then the default value is 7."""

    helices_view_order: List[int] = None
    """A list of the order in which the helix should be displayed in the main view of scadnano.
       
    This list must be a permutation containing each integer 0, 1, 2, ..., len(helices)-1 exactly once.

    Optional field. If not specified, it will be set to the identity permutation [0, ..., len(helices)-1].
    """

    @staticmethod
    def from_file(filename: str) -> DNADesign:
        """
        Loads a :any:`DNADesign` from the file with the given name.

        :param filename: name of the file with the design. Should be a JSON file ending in .dna
        :return: DNADesign described in the file
        """
        with open(filename) as f:
            json_str = f.read()
        json_map = json.loads(json_str)
        return DNADesign.from_json(json_map)

    @staticmethod
    def from_json(json_map: dict) -> DNADesign:
        version = json_map.get(version_key, initial_version)  # not sure what to do with this
        grid = json_map.get(grid_key, Grid.square)
        grid_is_none = grid == Grid.none

        if (major_tick_distance_key in json_map):
            major_tick_distance = json_map[major_tick_distance_key]
        elif not grid_is_none:
            if grid in [Grid.hex, Grid.honeycomb]:
                major_tick_distance = 7
            else:
                major_tick_distance = 8
        else:
            major_tick_distance = -1

        helices = []
        deserialized_helices_list = json_map[helices_key]
        num_helices = len(deserialized_helices_list)

        # create Helices
        idx = 0
        for helix_json in deserialized_helices_list:
            helix = Helix.from_json(helix_json)
            helix.set_idx(idx)
            if (grid_is_none and grid_position_key in helix_json):
                raise IllegalDNADesignError(
                    f'grid is none, but Helix {idx} has grid_position = {helix_json[grid_position_key]}')
            elif not grid_is_none and position3d_key in helix_json:
                raise IllegalDNADesignError(
                    'grid is not none, but Helix $idx has position = ${helix_json[constants.position3d_key]}')
            helices.append(helix)
            idx += 1

        # view order of helices
        helices_view_order = json_map.get(helices_view_order_key)
        if helices_view_order is not None:
            identity_permutation = list(range(num_helices))
            if len(helices_view_order) != num_helices:
                raise IllegalDNADesignError(f'length of helices ({num_helices}) does not match '
                                            f'length of helices_view_order ({len(helices_view_order)})')
            if sorted(helices_view_order) != identity_permutation:
                raise IllegalDNADesignError(f'helices_view_order = {helices_view_order} is not a permutation')

        # strands
        strands = []
        deserialized_strand_list = json_map[strands_key]
        for strand_json in deserialized_strand_list:
            strand = Strand.from_json(strand_json)
            strands.append(strand)

        return DNADesign(
            helices=helices,
            strands=strands,
            grid=grid,
            major_tick_distance=major_tick_distance,
            helices_view_order=helices_view_order,
        )

    @staticmethod
    def from_cadnano_v2(directory, filename) -> DNAOrigamiDesign:
        """ Creates a DNAOrigamiDesign from a cadnano v2 file.
        """
        file_path = os.path.join(directory, filename)
        f = open(file_path, 'r')
        cadnano_v2_design = json.load(f)
        f.close()

        num_bases = len(cadnano_v2_design['vstrands'][0]['scaf'])
        grid_type = Grid.square
        if num_bases % 21 == 0:
            grid_type = Grid.honeycomb
            raise NotImplementedError("Can't import honeycomb yet")

        min_row, min_col = None, None
        for cadnano_helix in cadnano_v2_design['vstrands']:
            col, row = cadnano_helix['col'], cadnano_helix['row']
            min_row = row if min_row is None else min_row
            min_col = col if min_col is None else min_col
            min_row = row if row < min_row else min_row
            min_col = col if col < min_col else min_col

        helices = []
        for cadnano_helix in cadnano_v2_design['vstrands']:
            col, row = cadnano_helix['col'], cadnano_helix['row']
            helix = Helix(max_offset=num_bases, grid_position=[col - min_col, row - min_row, 0])
            helix.set_idx(cadnano_helix['num'])
            helices.append(helix)

        sorted_helices = sorted(helices, key=lambda
            h: h.idx())  # Needed to show helices in the same order than cadnano in the grid view
        design = DNAOrigamiDesign(grid=grid_type, helices=sorted_helices, strands=[])
        design.set_helices_view_order([h.idx() for h in helices])

        return design

    def to_json_serializable(self, suppress_indent=True):
        dct = OrderedDict()
        dct[version_key] = current_version
        if self.grid != default_grid:
            dct[grid_key] = str(self.grid)[5:]  # remove prefix 'Grid.'
        if self.major_tick_distance >= 0 and (
                self.major_tick_distance != default_major_tick_distance(self.grid)):
            dct[major_tick_distance_key] = self.major_tick_distance

        dct[helices_key] = [helix.to_json_serializable(suppress_indent) for helix in self.helices]

        default_helices_view_order = list(range(0, len(self.helices)))
        if self.helices_view_order != default_helices_view_order:
            dct[helices_view_order_key] = _NoIndent(self.helices_view_order)

        dct[strands_key] = [strand.to_json_serializable(suppress_indent) for strand in self.strands]

        for helix in self.helices:
            helix_json = dct[helices_key][helix.idx()].value  # get past NoIndent surrounding helix
            # XXX: no need to check here because key was already deleted by Helix.to_json_serializable
            # max_offset still needs to be checked here since it requires global knowledge of Strands
            # if 0 == helix_json[min_offset_key]:
            #     del helix_json[min_offset_key]
            max_offset = max((ss.end for ss in helix._substrands), default=-1)
            if max_offset == helix_json[max_offset_key]:
                del helix_json[max_offset_key]

        return dct

    def _get_multiple_of_x_sup_closest_to_y(self, x: int, y: int) -> int:
        return y if y % x == 0 else y + (x - y % x)

    def _cadnano_v2_convert_honeycomb_coords(self, grid_position: Tuple[int, int, int]):
        """Converts scadnano honeycomb coords ("odd-r horizontal layout") to cadnano v2 convention. """
        raise NotImplementedError('The honeycomb lattice is not exportable yet')
        pass

    def _cadnano_v2_place_strand_segment(self, helix_dct, substrand: Substrand,
                                         strand_type: str = 'scaf') -> None:
        """Converts a strand region with no crossover to cadnano v2.
        """
        # Insertions and deletions
        for deletion in substrand.deletions:
            helix_dct['skip'][deletion] = -1
        for loop_where, loop_nb in substrand.insertions:
            helix_dct['loop'][loop_where] = loop_nb

        start, end, forward = substrand.start, substrand.end, substrand.forward
        strand_helix = helix_dct['num']

        for i_base in range(start, end):
            if forward:
                from_helix, from_base = strand_helix, i_base - 1
                to_helix, to_base = strand_helix, i_base + 1
            else:
                from_helix, from_base = strand_helix, i_base + 1
                to_helix, to_base = strand_helix, i_base - 1

            if i_base == start:
                if forward:
                    helix_dct[strand_type][i_base][2:] = [to_helix, to_base]
                else:
                    helix_dct[strand_type][i_base][:2] = [from_helix, from_base]
            elif i_base < end - 1:
                helix_dct[strand_type][i_base] = [from_helix, from_base, to_helix, to_base]
            else:
                if forward:
                    helix_dct[strand_type][i_base][:2] = [from_helix, from_base]
                else:
                    helix_dct[strand_type][i_base][2:] = [to_helix, to_base]
        return

    def _cadnano_v2_place_crossover(self, helix_from_dct, helix_to_dct,
                                    substrand_from: Substrand, substrand_to: Substrand,
                                    strand_type: str = 'scaf') -> None:
        """Converts a crossover to cadnano v2 format.
        Returns a conversion table from ids in the structure self.helices to helices ids
        as given by helix.idx().
        """

        helix_from = helix_from_dct['num']
        start_from, end_from, forward_from = substrand_from.start, substrand_from.end, substrand_from.forward

        helix_to = helix_to_dct['num']
        start_to, end_to = substrand_to.start, substrand_to.end

        if forward_from:
            helix_from_dct[strand_type][end_from - 1][2:] = [helix_to, end_to - 1]
            helix_to_dct[strand_type][end_to - 1][:2] = [helix_from, end_from - 1]
        else:
            helix_from_dct[strand_type][start_from][2:] = [helix_to, start_to]
            helix_to_dct[strand_type][start_to][:2] = [helix_from, start_from]

    def _cadnano_v2_color_of_stap(self, color, substrand) -> List[int]:
        base_id = substrand.start if substrand.forward else substrand.end - 1
        cadnano_color = color.to_cadnano_v2_int_hex()
        return [base_id, cadnano_color]

    def _cadnano_v2_place_strand(self, strand, dct, helices_ids_reverse) -> None:
        """Place a scadnano strand in cadnano v2.
        """
        strand_type = 'stap'
        if hasattr(strand, is_scaffold_key):
            strand_type = 'scaf'

        for i, substrand in enumerate(strand.substrands):
            which_helix_id = helices_ids_reverse[substrand.helix]
            which_helix = dct['vstrands'][which_helix_id]

            if strand_type == 'stap':
                which_helix['stap_colors'].append(self._cadnano_v2_color_of_stap(strand.color, substrand))

            self._cadnano_v2_place_strand_segment(which_helix, substrand, strand_type)

            if i != len(strand.substrands) - 1:
                next_substrand = strand.substrands[i + 1]
                next_helix_id = helices_ids_reverse[next_substrand.helix]
                next_helix = dct['vstrands'][next_helix_id]
                self._cadnano_v2_place_crossover(which_helix, next_helix,
                                                 substrand, next_substrand, strand_type)

    def _cadnano_v2_fill_blank(self, dct, num_bases) -> None:
        """Creates blank cadnanov2 helices in and initialized all their fields.
        """
        helices_ids_reverse = {}
        for i, helix in enumerate(self.helices):
            helix_dct = OrderedDict()
            helix_dct['num'] = helix.idx()

            if self.grid == Grid.square:
                helix_dct['row'] = helix.grid_position[1]
                helix_dct['col'] = helix.grid_position[0]

            if self.grid == Grid.honeycomb:
                helix_dct['row'], helix_dct['col'] = self._cadnano_v2_convert_honeycomb_coords(
                    helix.grid_position)

            helix_dct['scaf'] = []
            helix_dct['stap'] = []
            helix_dct['loop'] = []
            helix_dct['skip'] = []

            for _ in range(num_bases):
                helix_dct['scaf'].append([-1, -1, -1, -1])
                helix_dct['stap'].append([-1, -1, -1, -1])
                helix_dct['loop'].append(0)
                helix_dct['skip'].append(0)

            helix_dct['stap_colors'] = []
            helix_dct['scafLoop'] = []
            helix_dct['stapLoop'] = []

            helices_ids_reverse[helix_dct['num']] = i
            dct['vstrands'].append(helix_dct)
        return helices_ids_reverse

    def to_cadnano_v2(self):
        """Converts the design to the cadnano v2 format.
        Please see the spec `misc/cadnano-format-specs/v2.txt` for more info on that format.
        """
        dct = OrderedDict()
        dct['vstrands'] = []

        if self.__class__ != DNAOrigamiDesign:
            raise ValueError(
                'Please export DNAOrigamiDesign only as we need to know which strand is the scaffold.')

        '''Figuring out the type of grid.
        In cadnano v2, all helices have the same max offset 
        called `num_bases` and the type of grid is determined as follows:
            if num_bases % 32 == 0: then we are on grid square
            if num_bases % 21 == 0: then we are on grid honey
        '''
        num_bases = 0
        for helix in self.helices:
            num_bases = max(num_bases, helix.max_offset)

        if self.grid == Grid.square:
            num_bases = self._get_multiple_of_x_sup_closest_to_y(32, num_bases)
        elif self.grid == Grid.honeycomb:
            num_bases = self._get_multiple_of_x_sup_closest_to_y(21, num_bases)
        else:
            raise NotImplementedError('We can export to cadnano v2 `square` and `honeycomb` grids only.')

        '''Figuring out if helices numbers have good parity.
        In cadnano v2, only even helices have the scaffold go forward, only odd helices
        have the scaffold go backward.
        TODO: test for that case
        '''
        for strand in self.strands:
            if hasattr(strand, is_scaffold_key):
                for substrand in strand.substrands:
                    if type(substrand) == Loopout:
                        raise ValueError(
                            'We cannot handle designs with Loopouts as it is not a cadnano v2 concept')
                    if substrand.helix % 2 != int(not substrand.forward):
                        raise ValueError('We can only convert designs where even helices have the scaffold \
                                                  going forward and odd helices have the scaffold going backward see the spec v2.txt Note 4.')

        '''Filling the helices with blank.
        '''
        helices_ids_reverse = self._cadnano_v2_fill_blank(dct, num_bases)
        '''Putting the scaffold in place.
        '''

        for strand in self.strands:
            self._cadnano_v2_place_strand(strand, dct, helices_ids_reverse)

        return dct

    def __post_init__(self):
        if self.major_tick_distance < 0 or self.major_tick_distance is None:
            self.major_tick_distance = default_major_tick_distance(self.grid)

        # doing this first matters because most of DNADesign assumes helices has been set
        if self.helices is None:
            if len(self.strands) > 0:
                max_helix_idx = max(ss.helix for strand in self.strands for ss in strand.bound_substrands())
                self.helices = list(Helix() for _ in range(max_helix_idx + 1))
            else:
                self.helices = []

        # XXX: exact order of these calls is important
        self._set_helices_idxs()
        self._set_helices_grid_and_svg_positions()
        self._build_substrands_on_helix_lists()
        self._set_helices_min_max_offsets(update=False)
        self._check_legal_design()

        self._set_and_check_helices_view_order()

    def _set_and_check_helices_view_order(self):
        identity = list(range(0, len(self.helices)))
        if self.helices_view_order is None:
            self.helices_view_order = identity
        self._check_helices_view_order_is_bijection()

    def set_helices_view_order(self, helices_view_order: List[int]):
        self.helices_view_order = helices_view_order
        self._check_helices_view_order_is_bijection()

    def _check_helices_view_order_is_bijection(self):
        identity = list(range(0, len(self.helices)))
        if not (sorted(self.helices_view_order) == identity):
            raise IllegalDNADesignError(
                f"The specified helices view order: {self.helices_view_order}\n "
                f"is not a bijection from [0,{len(self.helices) - 1}] to [0,{len(self.helices) - 1}].")

    def _set_helices_idxs(self):
        for idx, helix in enumerate(self.helices):
            helix.set_idx(idx)

    def _set_helices_grid_and_svg_positions(self):
        for idx, helix in enumerate(self.helices):
            if helix.grid_position is None:
                helix.grid_position = helix.default_grid_position()
            if helix.svg_position is None:
                helix.svg_position = helix.default_svg_position()

    def _set_helices_min_max_offsets(self, update: bool):
        """update = whether to overwrite existing Helix.max_offset and Helix.min_offset.
        Don't do this when DNADesign is first created, but do it later when updating."""
        for helix in self.helices:

            if update or helix.max_offset is None:
                max_offset = None if len(helix._substrands) == 0 else helix._substrands[0].end
                for substrand in helix._substrands:
                    max_offset = max(max_offset, substrand.end)
                helix.max_offset = max_offset

            if update or helix.min_offset is None:
                min_offset = None if len(helix._substrands) == 0 else helix._substrands[0].start
                for substrand in helix._substrands:
                    min_offset = min(min_offset, substrand.start)
                if min_offset > 0: min_offset = 0
                helix.min_offset = min_offset

    def set_default_idt(self, use_default_idt):
        """If ``True``, sets :py:data:`Strand.use_default_idt` to ``True`` for every :any:`Strand` in this
        :any:`DNADesign` and calls :py:meth:`Strand.set_default_idt` on each of them to assign a
        default idt field.

        If ``False``, removes IDT field from each :any:`Strand`."""
        for strand in self.strands:
            strand.set_default_idt(use_default_idt)

    def strands_starting_on_helix(self, helix: int) -> List[Strand]:
        """Return list of :any:`Strand`'s that begin (have their 5' end)
        on the :any:`Helix` with index `helix`."""
        return [strand for strand in self.strands if strand.substrands[0].helix == helix]

    def strands_ending_on_helix(self, helix: int) -> List[Strand]:
        """Return list of :any:`Strand`'s that finish (have their 3' end)
        on the :any:`Helix` with index `helix`."""
        return [strand for strand in self.strands if strand.substrands[-1].helix == helix]

    def _check_legal_design(self):
        # self._check_helix_indices()
        self._check_helix_offsets()
        self._check_strands_reference_helices_legally()
        self._check_loopouts_not_consecutive_or_singletons_or_zero_length()
        self._check_strands_overlap_legally()
        self._check_grid_honeycomb_positions_legal()

    def _check_grid_honeycomb_positions_legal(self):
        # ensures grid positions are legal if honeycomb lattice is used
        if self.grid == Grid.honeycomb:
            for helix in self.helices:
                x = helix.grid_position[0]
                y = helix.grid_position[1]

                # following is for odd-q system: https://www.redblobgames.com/grids/hexagons/
                if x % 2 == 1 and y % 2 == 0:
                    raise IllegalDNADesignError('honeycomb lattice disallows grid positions of first two '
                                                'coordinates (x,y,_) if x is odd and y is even, '
                                                f'but helix {helix.idx()} has grid position '
                                                f'{helix.grid_position}')

                # following is for even-q system: https://www.redblobgames.com/grids/hexagons/
                # if x % 2 == 1 and y % 2 == 1:
                #     raise IllegalDNADesignError('honeycomb lattice disallows grid positions of first two '
                #                                 'coordinates (x,y,_) if both x and y are odd, '
                #                                 f'but helix {helix.idx()} has grid position '
                #                                 f'{helix.grid_position}')

                # following is for odd-r system: https://www.redblobgames.com/grids/hexagons/
                # if x % 3 == 0 and y % 2 == 0:
                #     raise IllegalDNADesignError('honeycomb lattice disallows grid positions of first two '
                #                                 'coordinates (x,y,_) with y even and x a multiple of 3, '
                #                                 f'but helix {helix.idx()} has grid position '
                #                                 f'{helix.grid_position}')
                # if x % 3 == 1 and y % 2 == 1:
                #     raise IllegalDNADesignError('honeycomb lattice disallows grid positions of first two '
                #                                 'coordinates (x,y) with y odd and x = 1 + a multiple of 3, '
                #                                 f'but helix {helix.idx()} has grid position '
                #                                 f'{helix.grid_position}')

    # TODO: come up with reasonable default behavior when no strands are on helix and max_offset not given
    def _check_helix_offsets(self):
        for helix in self.helices:
            if helix.min_offset >= helix.max_offset:
                err_msg = f'for helix {helix.idx()}, ' \
                          f'helix.min_offset = {helix.min_offset} must be strictly less than ' \
                          f'helix.max_offset = {helix.max_offset}'
                raise IllegalDNADesignError(err_msg)

    # def _check_helix_indices(self):
    #     # ensure if there are H helices, the list of sorted indices is 0,1,...,H-1
    #     indices_helices = sorted([(helix.idx, helix) for helix in self.helices],
    #                              key=lambda x: x[0])
    #     for (correct_idx, (helix_idx, helix)) in enumerate(indices_helices):
    #         if correct_idx != helix_idx:
    #             if correct_idx < helix_idx:
    #                 err_msg = f"missing Helix with helix {correct_idx}"
    #             else:
    #                 err_msg = f"duplicate Helices with helix {helix_idx}"
    #             raise IllegalDNADesignError(err_msg)

    def _check_strands_overlap_legally(self, substrand_to_check: Substrand = None):
        """If `substrand_to_check` is None, check all.
        Otherwise only check pairs where one is substrand_to_check."""

        def err_msg(substrand1, substrand2, h_idx):
            return f"two substrands overlap on helix {h_idx}: " \
                   f"\n{substrand1}\n  and\n{substrand2}\n  but have the same direction"

        # ensure that if two strands overlap on the same helix,
        # they point in opposite directions
        for helix_idx, substrands in enumerate(helix._substrands for helix in self.helices):
            if substrand_to_check is not None and substrand_to_check.helix != helix_idx:
                # TODO: if necessary, we can be more efficient by only checking this one substrand
                continue

            if len(substrands) == 0:
                continue

            # check all consecutive substrands on the same helix, sorted by start/end indices
            offsets_data = []
            for substrand in substrands:
                offsets_data.append((substrand.start, True, substrand))
                offsets_data.append((substrand.end, False, substrand))
            offsets_data.sort(key=lambda offset_data: offset_data[0])

            current_substrands: List[Substrand] = []
            for offset, is_start, substrand in offsets_data:
                if is_start:
                    if len(current_substrands) >= 2:
                        if offset >= current_substrands[1].end:
                            del current_substrands[1]
                    if len(current_substrands) >= 1:
                        if offset >= current_substrands[0].end:
                            del current_substrands[0]
                    current_substrands.append(substrand)
                    if len(current_substrands) > 2:
                        ss0, ss1, ss2 = current_substrands[0:3]
                        for s_first, s_second in [(ss0, ss1), (ss1, ss2), (ss0, ss2)]:
                            if s_first.forward == s_second.forward:
                                raise IllegalDNADesignError(err_msg(s_first, s_second, helix_idx))
                        raise AssertionError(
                            f"since current_substrands = {current_substrands} has at least three substrands, "
                            f"I expected to find a pair of illegally overlapping substrands")
                    elif len(current_substrands) == 2:
                        s_first, s_second = current_substrands
                        if s_first.forward == s_second.forward:
                            raise IllegalDNADesignError(err_msg(s_first, s_second, helix_idx))

    def _check_loopouts_not_consecutive_or_singletons_or_zero_length(self):
        for strand in self.strands:
            DNADesign._check_loopout_not_singleton(strand)
            DNADesign._check_two_consecutive_loopouts(strand)
            DNADesign._check_loopouts_length(strand)

    @staticmethod
    def _check_two_consecutive_loopouts(strand):
        for ss1, ss2 in _pairwise(strand.substrands):
            if ss1.is_loopout() and ss2.is_loopout():
                raise StrandError(strand, 'cannot have two consecutive Loopouts in a strand')

    @staticmethod
    def _check_loopout_not_singleton(strand):
        if len(strand.substrands) == 1:
            if strand.first_substrand().is_loopout():
                raise StrandError(strand, 'strand cannot have a single Loopout as its only substrand')

    @staticmethod
    def _check_loopouts_length(strand):
        for loopout in strand.substrands:
            if loopout.is_loopout() and loopout.length <= 0:
                raise StrandError(strand, f'loopout length must be positive but is {loopout.length}')

    def _check_strands_reference_helices_legally(self):
        # ensure each strand refers to an existing helix
        for strand in self.strands:
            self._check_strand_references_legal_helices(strand)
            self._check_strand_has_legal_offsets_in_helices(strand)

    def _check_strand_has_legal_offsets_in_helices(self, strand: Strand):
        for substrand in strand.substrands:
            if substrand.is_substrand():
                helix = self.helices[substrand.helix]
                if substrand.start < helix.min_offset:
                    err_msg = f"substrand {substrand} has start offset {substrand.start}, " \
                              f"beyond the end of " \
                              f"Helix {substrand.helix} that has min_offset = {helix.min_offset}"
                    raise StrandError(substrand._parent_strand, err_msg)
                if substrand.end > helix.max_offset:
                    err_msg = f"substrand {substrand} has end offset {substrand.end}, " \
                              f"beyond the end of " \
                              f"Helix {substrand.helix} that has max_offset = {helix.max_offset}"
                    err = StrandError(substrand._parent_strand, err_msg)
                    raise err

    def _check_strand_references_legal_helices(self, strand: Strand):
        for substrand in strand.substrands:
            if substrand.is_substrand() and not (0 <= substrand.helix < len(self.helices)):
                err_msg = f"substrand {substrand} refers to nonexistent Helix index {substrand.helix}"
                raise StrandError(substrand._parent_strand, err_msg)

        # ensure helix_idx's are never negative twice in a row
        for ss1, ss2 in _pairwise(strand.substrands):
            if ss1.is_loopout() and ss2.is_loopout():
                err_msg = f"Loopouts {ss1} and {ss2} are consecutive on strand {strand}. " \
                          f"At least one of any consecutive pair must be a Substrand, not a Loopout."
                raise StrandError(strand, err_msg)

    def substrand_at(self, helix: int, offset: int, forward: bool):
        """
        Return :any:`Substrand` that overlaps `offset` on helix with idx `helix` and has
        :py:data:`Substrand.forward` = ``True``, or ``None`` if there is no such :any:`Substrand`.

        :param helix: TODO
        :param offset: TODO
        :param forward: TODO
        :return: TODO
        """
        for substrand in self.substrands_at(helix, offset):
            if substrand.forward == forward:
                return substrand
        return None

    def substrands_at(self, helix: int, offset: int) -> List[Substrand]:
        """Return list of :any:`Substrand`'s that overlap `offset` on helix with idx `helix`.

        If constructed properly, this list should have 0, 1, or 2 elements."""
        substrands_on_helix = self.helices[helix]._substrands
        # TODO: replace this with a faster algorithm using binary search
        substrands_on_helix = [substrand for substrand in substrands_on_helix if
                               substrand.contains_offset(offset)]
        if len(substrands_on_helix) not in [0, 1, 2]:
            raise AssertionError(f'There should be at most 2 substrands on helix {helix}, '
                                 f'but there are {len(substrands_on_helix)}:\n{substrands_on_helix}')
        return substrands_on_helix

    # TODO: add_strand and insert_substrand should check for existing deletions/insertion parallel strands
    def add_strand(self, strand: Strand):
        """Add `strand` to this design."""
        self._check_strand_references_legal_helices(strand)
        self.strands.append(strand)
        for substrand in strand.substrands:
            if substrand.is_substrand():
                self.helices[substrand.helix]._substrands.append(substrand)

    def remove_strand(self, strand: Strand):
        """Remove `strand` from this design."""
        self.strands.remove(strand)
        for substrand in strand.substrands:
            if substrand.is_substrand():
                self.helices[substrand.helix]._substrands.remove(substrand)

    def insert_substrand(self, strand: Strand, order: int, substrand: Union[Substrand, Loopout]):
        """Insert `substrand` into `strand` at index given by `order`. Uses same indexing as Python lists,
        e.g., ``strand.insert_substrand(ss, 0)`` inserts ``ss`` as the new first :any:`Substrand`."""
        assert strand in self.strands
        strand._insert_substrand(order, substrand)
        self._check_strand_references_legal_helices(strand)
        self._check_loopouts_not_consecutive_or_singletons_or_zero_length()
        if substrand.is_substrand():
            self.helices[substrand.helix]._substrands.append(substrand)

    def remove_substrand(self, strand: Strand, substrand: Union[Substrand, Loopout]):
        """Remove `substrand` from `strand`."""
        assert strand in self.strands
        strand._remove_substrand(substrand)
        if substrand.is_substrand():
            self.helices[substrand.helix]._substrands.remove(substrand)

    def _build_substrands_on_helix_lists(self):
        for helix in self.helices:
            helix._substrands = []
        for strand in self.strands:
            for substrand in strand.substrands:
                if substrand.is_substrand():
                    if substrand.helix < len(self.helices):
                        self.helices[substrand.helix]._substrands.append(substrand)
                    else:
                        msg = f"substrand's helix is {substrand.helix} but largest helix is {len(self.helices) - 1}"
                        raise StrandError(strand=strand, the_cause=msg)

    def to_json(self, suppress_indent=True) -> str:
        """Return string representing this DNADesign, suitable for reading by scadnano if written to
        a JSON file ending in extension .dna"""
        return _json_encode(self, suppress_indent)

    # TODO: create version of add_deletion and add_insertion that simply changes the major tick distance
    #  on the helix at that position, as well as updating the end offset of the substrand (and subsequent
    #  substrands on the same helix)

    def add_deletion(self, helix: int, offset: int):
        """Adds a deletion to every :class:`scadnano.Strand` at the given helix and base offset."""
        substrands = self.substrands_at(helix, offset)
        if len(substrands) == 0:
            raise IllegalDNADesignError(f"no substrands are at helix {helix} offset {offset}")
        for substrand in substrands:
            if substrand.contains_offset(offset):
                substrand.deletions.append(offset)

    def add_insertion(self, helix: int, offset: int, length: int):
        """Adds an insertion with the given length to every :class:`scadnano.Strand`
        at the given helix and base offset, with the given length."""
        substrands = self.substrands_at(helix, offset)
        if len(substrands) == 0:
            raise IllegalDNADesignError(f"no substrands are at helix {helix} offset {offset}")
        for substrand in substrands:
            if substrand.contains_offset(offset):
                substrand.insertions.append((offset, length))

    def set_start(self, substrand: Substrand, start: int):
        """Sets ``substrand.start`` to `start`."""
        assert substrand in (ss for strand in self.strands for ss in strand.substrands)
        substrand.set_start(start)
        self._check_strands_overlap_legally(substrand)

    def set_end(self, substrand: Substrand, end: int):
        """Sets ``substrand.end`` to `end`."""
        assert substrand in (ss for strand in self.strands for ss in strand.substrands)
        substrand.set_end(end)
        self._check_strands_overlap_legally(substrand)

    def move_strand_offsets(self, delta: int):
        """Moves all strands backward (if `delta` < 0) or forward (if `delta` > 0) by `delta`."""
        for strand in self.strands:
            for substrand in strand.substrands:
                substrand.start += delta
                substrand.end += delta
        self._check_strands_overlap_legally()

    def move_strands_on_helices(self, delta: int):
        """Moves all strands up (if `delta` < 0) or down (if `delta` > 0) by the number of helices given by
        `delta`."""
        for strand in self.strands:
            for substrand in strand.substrands:
                substrand.helix += delta
        self._check_strands_reference_helices_legally()

    def assign_dna(self, strand: Strand, sequence: str):
        """
        Assigns `sequence` as DNA sequence of `strand`.

        If any :class:`scadnano.Strand` is bound to `strand`,
        it is assigned the reverse Watson-Crick complement of the relevant portion,
        and any remaining portions of the other strand that have not already been assigned a DNA sequence
        are assigned to be the symbol :py:data:`DNA_base_wildcard`.

        Before assigning, `sequence` is first forced to be the same length as `strand`
        as follows:
        If `sequence` is longer, it is truncated.
        If `sequence` is shorter, it is padded with :py:data:`DNA_base_wildcard`'s.

        All whitespace in `sequence` is removed, and lowercase bases
        'a', 'c', 'g', 't' are converted to uppercase.
        """
        padded_sequence = _pad_and_remove_whitespace(sequence, strand)
        if strand is None:
            raise IllegalDNADesignError('strand cannot be None to assign DNA to it')
        if strand not in self.strands:
            raise StrandError(strand, 'strand is not in the given DNADesign')

        if strand.dna_sequence is None:
            merged_sequence = padded_sequence
        else:
            try:
                merged_sequence = _string_merge_wildcard(strand.dna_sequence, padded_sequence,
                                                         DNA_base_wildcard)
            except ValueError:
                first_ss = strand.first_substrand()
                msg = f'strand starting at helix {first_ss.helix}, offset {first_ss.offset_5p()} has ' \
                      f'length ' \
                      f'{strand.dna_length()} and already has a DNA sequence assignment of length ' \
                      f'{len(strand.dna_sequence)}, which is \n' \
                      f'{strand.dna_sequence}, ' \
                      f'but you tried to assign a different sequence of length {len(padded_sequence)} to ' \
                      f'it, which is\n{padded_sequence}.'
                raise IllegalDNADesignError(msg)

        strand.set_dna_sequence(merged_sequence)

        for other_strand in self.strands:
            # note that possibly strand==other_strand; it might bind to itself at some point and we want to
            # allow a partial assignment to one substrand to automatically assign the complement to the
            # bound substrand.
            # However, if there are no wildcards in the assigned sequence we can safely skip strand.
            if strand == other_strand and DNA_base_wildcard not in strand.dna_sequence:
                continue
            if other_strand.overlaps(strand):
                # we do this even if other_strand has a complete DNA sequence,
                # because we get complementarity checks this way
                other_strand.assign_dna_complement_from(strand)

    def to_idt_bulk_input_format(self, delimiter: str = ',', warn_duplicate_name: bool = False,
                                 warn_on_non_idt_strands: bool = False) -> str:
        """Return string that is written to the file in the method
        :py:meth:`DNADesign.write_idt_bulk_input_file`.

        `delimiter` is the symbol to delimit the four IDT fields name,sequence,scale,purification.

        `warn_duplicate_name` if ``True`` prints a warning when two different :any:`Strand`'s have the same
        :py:attr:`IDTField.name` and the same :any:`Strand.dna_sequence`. An :any:`IllegalDNADesignError` is
        raised (regardless of the value of this parameter)
        if two different :any:`Strand`'s have the same name but different sequences, IDT scales, or IDT
        purifications.

        `warn_on_non_idt_strands` specifies whether to print a warning for strands that lack the field
        :any:`Strand.idt`. Such strands will not be part of the output.
        """
        added_strands = self._idt_strands(warn_duplicate_name, warn_on_non_idt_strands)

        idt_lines = [
            delimiter.join([strand.idt.name, strand.dna_sequence, strand.idt.scale, strand.idt.purification])
            for strand in added_strands.values()]

        idt_string = '\n'.join(idt_lines)
        return idt_string

    def _idt_strands(self, warn_duplicate_name, warn_on_non_idt_strands) -> Dict[str, Strand]:
        added_strands: Dict[str, Strand] = {}  # dict: name -> strand
        for strand in self.strands:
            if strand.idt is not None:
                name = strand.idt.name
                if name in added_strands:
                    existing_strand = added_strands[name]
                    assert existing_strand.idt.name == name
                    ss = strand.first_substrand()
                    existing_ss = existing_strand.first_substrand()
                    if strand.dna_sequence != existing_strand.dna_sequence:
                        raise IllegalDNADesignError(
                            f'two strands with same IDT name {name} but different sequences:\n'
                            f'  strand 1: helix {ss.helix}, 5\' end at offset {ss.offset_5p()}, '
                            f'sequence: {strand.dna_sequence}\n'
                            f'  strand 2: helix {existing_ss.helix}, 5\' end at offset '
                            f'{existing_ss.offset_5p()}, '
                            f'sequence: {existing_strand.dna_sequence}\n')
                    elif strand.idt.scale != existing_strand.idt.scale:
                        raise IllegalDNADesignError(
                            f'two strands with same IDT name {name} but different IDT scales:\n'
                            f'  strand 1: helix {ss.helix}, 5\' end at offset {ss.offset_5p()}, '
                            f'scale: {strand.idt.scale}\n'
                            f'  strand 2: helix {existing_ss.helix}, 5\' end at offset '
                            f'{existing_ss.offset_5p()}, '
                            f'scale: {existing_strand.idt.scale}\n')
                    elif strand.idt.purification != existing_strand.idt.purification:
                        raise IllegalDNADesignError(
                            f'two strands with same IDT name {name} but different purifications:\n'
                            f'  strand 1: helix {ss.helix}, 5\' end at offset {ss.offset_5p()}, '
                            f'purification: {strand.idt.purification}\n'
                            f'  strand 2: helix {existing_ss.helix}, 5\' end at offset '
                            f'{existing_ss.offset_5p()}, '
                            f'purification: {existing_strand.idt.purification}\n')
                    elif warn_duplicate_name:
                        print(
                            f'WARNING: two strands with same IDT name {name}:\n'
                            f'  strand 1: helix {ss.helix}, 5\' end at offset {ss.offset_5p()}\n'
                            f'  strand 2: helix {existing_ss.helix}, 5\' end at offset '
                            f'{existing_ss.offset_5p()}\n')
                added_strands[name] = strand
            elif warn_on_non_idt_strands:
                print(f"WARNING: strand with 5' end on helix {strand.first_substrand().helix} "
                      f"does not have a field idt, so will not be part of IDT output.")
        return added_strands

    def write_idt_bulk_input_file(self, directory: str = '.', filename=None, delimiter: str = ',',
                                  warn_duplicate_name: bool = False, warn_on_non_idt_strands=False):
        """Write ``.idt`` text file encoding the strands of this :any:`DNADesign` with the field
        :any:`Strand.idt`, suitable for pasting into the "Bulk Input" field of IDT
        (Integrated DNA Technologies, Coralville, IA, https://www.idtdna.com/),
        with the output file having the same name as the running script but with ``.py`` changed to ``.idt``,
        unless `filename` is explicitly specified.
        For instance, if the script is named ``my_origami.py``,
        then the sequences will be written to ``my_origami.idt``.

        `directory` specifies a directory in which to place the file, either absolute or relative to
        the current working directory. Default is the current working directory.

        `delimiter` is the symbol to delimit the four IDT fields name,sequence,scale,purification.

        `warn_duplicate_name` if ``True`` prints a warning when two different :any:`Strand`'s have the same
        :py:attr:`IDTField.name` and the same :any:`Strand.dna_sequence`. An :any:`IllegalDNADesignError` is
        raised (regardless of the value of this parameter)
        if two different :any:`Strand`'s have the same name but different sequences, IDT scales, or IDT
        purifications.

        `warn_on_non_idt_strands` specifies whether to print a warning for strands that lack the field
        :any:`Strand.idt`. Such strands will not be output into the file.

        The string written is that returned by :meth:`DNADesign.to_idt_bulk_input_format`.
        """
        contents = self.to_idt_bulk_input_format(delimiter, warn_duplicate_name, warn_on_non_idt_strands)
        _write_file_same_name_as_running_python_script(contents, 'idt', directory, filename)

    def write_idt_plate_excel_file(self, directory: str = '.', filename=None,
                                   warn_duplicate_name: bool = False, warn_on_non_idt_strands=False,
                                   use_default_plates=False, warn_using_default_plates=True,
                                   plate_type: PlateType = PlateType.wells96):
        """Write ``.xls`` (Microsoft Excel) file encoding the strands of this :any:`DNADesign` with the field
        :py:data:`Strand.idt`, suitable for uploading to IDT
        (Integrated DNA Technologies, Coralville, IA, https://www.idtdna.com/)
        to describe a 96-well or 384-well plate
        (https://www.idtdna.com/site/order/plate/index/dna/),
        with the output file having the same name as the running script but with ``.py`` changed to ``.xls``,
        unless `filename` is explicitly specified.
        For instance, if the script is named ``my_origami.py``,
        then the sequences will be written to ``my_origami.xls``.

        `directory` specifies a directory in which to place the file, either absolute or relative to
        the current working directory. Default is the current working directory.

        `warn_duplicate_name` if ``True`` prints a warning when two different :any:`Strand`'s have the same
        :py:attr:`IDTField.name` and the same :any:`Strand.dna_sequence`. An :any:`IllegalDNADesignError` is
        raised (regardless of the value of this parameter)
        if two different :any:`Strand`'s have the same name but different sequences, IDT scales, or IDT
        purifications.

        `warn_on_non_idt_strands` specifies whether to print a warning for strands that lack the field
        :any:`Strand.idt`. Such strands will not be output into the file.

        `plate_type` is a :any:`PlateType` specifying whether to use a 96-well plate or a 384-well plate
        if the `use_default_plates` parameter is ``True``.
        Ignored if `use_default_plates` is ``False``, because in that case the wells are explicitly set
        by the user, who is free to use coordinates for either plate type.
        """

        idt_strands = list(self._idt_strands(warn_duplicate_name, warn_on_non_idt_strands).values())

        if not use_default_plates:
            self._write_plates_assuming_explicit_in_each_strand(directory, filename, idt_strands)
        else:
            if warn_using_default_plates:
                print("WARNING: ignoring plate data in each strand and using default sequential assignment "
                      "of plates and wells")
            self._write_plates_default(directory, filename, idt_strands, plate_type=plate_type)

    def _write_plates_assuming_explicit_in_each_strand(self, directory: str, filename: str,
                                                       idt_strands: List[Strand]):
        plates = list({strand.idt.plate for strand in idt_strands if strand.idt is not None if
                       strand.idt.plate is not None})
        if len(plates) == 0:
            raise ValueError('Cannot write a a plate file since no plate data exists in any Strands '
                             'in the design.\n'
                             'Set the option use_default_plates=True in '
                             "DNADesign.write_idt_plate_excel_file\nif you don't want to enter plate "
                             'and well positions for each Strand you wish to write to the Excel file.')
        plates.sort()
        filename_plate, workbook = self._setup_excel_file(directory, filename)
        for plate in plates:
            worksheet = self._add_new_excel_plate_sheet(plate, workbook)

            strands_in_plate = [strand for strand in idt_strands if
                                strand.idt is not None and strand.idt.plate == plate]

            strands_in_plate.sort(key=lambda s: (int(s.idt.well[1:]), s.idt.well[0]))

            for row, strand in enumerate(strands_in_plate):
                worksheet.write(row + 1, 0, strand.idt.well)
                worksheet.write(row + 1, 1, strand.idt.name)
                worksheet.write(row + 1, 2, strand.dna_sequence)

            workbook.save(filename_plate)

    def _add_new_excel_plate_sheet(self, plate_name: str, workbook: xlwt.Workbook) -> xlwt.Worksheet:
        worksheet = workbook.add_sheet(plate_name)
        worksheet.write(0, 0, 'Well Position')
        worksheet.write(0, 1, 'Name')
        worksheet.write(0, 2, 'Sequence')
        return worksheet

    def _setup_excel_file(self, directory, filename):
        plate_extension = f'xls'
        if filename is None:
            filename_plate = _get_filename_same_name_as_running_python_script(
                directory, plate_extension, filename)
        else:
            filename_plate = _create_directory_and_set_filename(directory, filename)
        workbook = xlwt.Workbook()
        return filename_plate, workbook

    def _write_plates_default(self, directory: str, filename: str, idt_strands: List[Strand],
                              plate_type: PlateType = PlateType.wells96):
        plate_coord = _PlateCoordinate(plate_type=plate_type)
        plate = 1
        excel_row = 1
        filename_plate, workbook = self._setup_excel_file(directory, filename)
        worksheet = self._add_new_excel_plate_sheet(f'plate{plate}', workbook)

        for strand in idt_strands:
            well = plate_coord.well()
            worksheet.write(excel_row, 0, well)
            worksheet.write(excel_row, 1, strand.idt.name)
            worksheet.write(excel_row, 2, strand.dna_sequence)
            plate_coord.increment()
            if plate != plate_coord.plate():
                workbook.save(filename_plate)
                plate = plate_coord.plate()
                worksheet = self._add_new_excel_plate_sheet(f'plate{plate}', workbook)
                excel_row = 1
            else:
                excel_row += 1

        workbook.save(filename_plate)

    def write_scadnano_file(self, directory: str = '.', filename=None):
        """Write ``.dna`` file representing this :any:`DNADesign`, suitable for reading by scadnano,
        with the output file having the same name as the running script but with ``.py`` changed to ``.dna``,
        unless `filename` is explicitly specified.
        For instance, if the script is named ``my_origami.py``,
        then the design will be written to ``my_origami.dna``.

        `directory` specifies a directory in which to place the file, either absolute or relative to
        the current working directory. Default is the current working directory.

        The string written is that returned by :meth:`DNADesign.to_json`.
        """
        contents = self.to_json()
        _write_file_same_name_as_running_python_script(contents, 'dna', directory, filename)

    def export_cadnano_v2(self, directory: str = '.', filename=None):
        """Write ``.json`` file representing this :any:`DNADesign`, suitable for reading by cadnano v2,
        with the output file having the same name as the running script but with ``.py`` changed to ``.json``,
        unless `filename` is explicitly specified.
        For instance, if the script is named ``my_origami.py``,
        then the design will be written to ``my_origami.json``.

        `directory` specifies a directory in which to place the file, either absolute or relative to
        the current working directory. Default is the current working directory.

        The string written is that returned by :meth:`DNADesign.to_cadnano_v2`.
        """
        content_serializable = OrderedDict({})
        content_serializable['name'] = _get_filename_same_name_as_running_python_script(directory, 'json',
                                                                                        filename)
        content_serializable_final = self.to_cadnano_v2()
        content_serializable.update(content_serializable_final)

        encoder = _SuppressableIndentEncoder
        contents = json.dumps(content_serializable, cls=encoder, indent=2)

        _write_file_same_name_as_running_python_script(contents, 'json', directory, filename)

    def add_nick(self, helix: int, offset: int, forward: bool):
        """Add nick to :any:`Substrand` on :any:`Helix` with index `helix`,
        in direction given by `forward`, at offset `offset`. The two :any:`Substrand`'s created by this nick
        will have 5'/3' ends at offsets `offset` and `offset-1`.

        For example, if there is a :any:`Substrand` with
        :py:data:`Substrand.helix` = ``0``,
        :py:data:`Substrand.forward` = ``True``,
        :py:data:`Substrand.start` = ``0``,
        :py:data:`Substrand.end` = ``10``,
        then calling ``add_nick(helix=0, offset=5, forward=True)`` will split it into two :any:`Substrand`'s,
        with one substrands having the fields
        :py:data:`Substrand.helix` = ``0``,
        :py:data:`Substrand.forward` = ``True``,
        :py:data:`Substrand.start` = ``0``,
        :py:data:`Substrand.end` = ``5``,
        (recall that :py:data:`Substrand.end` is exclusive, meaning that the largest offset on this
        substrand is 4 = ``offset-1``)
        and the other substrand having the fields
        :py:data:`Substrand.helix` = ``0``,
        :py:data:`Substrand.forward` = ``True``,
        :py:data:`Substrand.start` = ``5``,
        :py:data:`Substrand.end` = ``10``.
        """
        for substrand_to_remove in self.substrands_at(helix, offset):
            if substrand_to_remove.forward == forward:
                break
        else:
            raise IllegalDNADesignError(f'no substrand at helix {helix} in direction '
                                        f'{"forward" if forward else "reverse"} at offset {offset}')
        strand = substrand_to_remove.strand()
        substrands = strand.substrands
        order = substrands.index(substrand_to_remove)
        substrands_before = substrands[:order]
        substrands_after = substrands[order + 1:]
        substrand_left = Substrand(helix, forward, substrand_to_remove.start, offset)
        substrand_right = Substrand(helix, forward, offset, substrand_to_remove.end)

        if substrand_to_remove.forward:
            substrand_to_add_before = substrand_left
            substrand_to_add_after = substrand_right
        else:
            substrand_to_add_before = substrand_right
            substrand_to_add_after = substrand_left

        if strand.dna_sequence:
            dna_sequence_before = ''.join(ss.dna_sequence() for ss in substrands_before)
            dna_sequence_after = ''.join(ss.dna_sequence() for ss in substrands_after)
            dna_sequence_on_substrand_left = substrand_to_remove.dna_sequence_in(substrand_to_remove.start,
                                                                                 offset - 1)
            dna_sequence_on_substrand_right = substrand_to_remove.dna_sequence_in(offset,
                                                                                  substrand_to_remove.end - 1)
            if substrand_to_remove.forward:
                dna_sequence_on_substrand_before = dna_sequence_on_substrand_left
                dna_sequence_on_substrand_after = dna_sequence_on_substrand_right
            else:
                dna_sequence_on_substrand_before = dna_sequence_on_substrand_right
                dna_sequence_on_substrand_after = dna_sequence_on_substrand_left
            dna_sequence_before_whole = dna_sequence_before + dna_sequence_on_substrand_before
            dna_sequence_after_whole = dna_sequence_on_substrand_after + dna_sequence_after
        else:
            dna_sequence_before_whole = None
            dna_sequence_after_whole = None

        self.strands.remove(strand)

        idt_present = strand.idt is not None
        strand_before = Strand(substrands=substrands_before + [substrand_to_add_before],
                               dna_sequence=dna_sequence_before_whole,
                               color=strand.color, idt=strand.idt if idt_present else None)

        strand_after = Strand(substrands=[substrand_to_add_after] + substrands_after,
                              dna_sequence=dna_sequence_after_whole,
                              automatically_assign_color=True, use_default_idt=idt_present)

        self.helices[helix]._substrands.remove(substrand_to_remove)
        self.helices[helix]._substrands.extend([substrand_to_add_before, substrand_to_add_after])

        self.strands.extend([strand_before, strand_after])

    def add_half_crossover(self, helix1: int, helix2: int, offset1: int, forward1: bool,
                           offset2: int = None, forward2: bool = None):
        """
        Add a half crossover from helix `helix1` at offset `offset1` to `helix2`, on the strand
        with :py:data:`Strand.forward` = `forward`.

        Unlike :py:meth:`DNADesign.add_full_crossover`, which automatically adds a nick between the two
        half-crossovers, to call this method, there must *already* be nicks adjacent to the given
        offsets on the given helices. (either on the left or right side)

        :param helix1: index of one helix of half crossover
        :param helix2: index of other helix of half crossover
        :param offset1: offset on `helix1` at which to add half crossover
        :param forward1: direction of :any:`Strand` on `helix1` to which to add half crossover
        :param offset2: offset on `helix2` at which to add half crossover.
            If not specified, defaults to `offset1`
        :param forward2: direction of :any:`Strand` on `helix2` to which to add half crossover.
            If not specified, defaults to the negation of `forward1`

        """
        if offset2 is None:
            offset2 = offset1
        if forward2 is None:
            forward2 = not forward1
        ss1 = self.substrand_at(helix1, offset1, forward1)
        ss2 = self.substrand_at(helix2, offset2, forward2)
        if ss1 is None:
            raise IllegalDNADesignError(f"Cannot add half crossover at (helix={helix1}, offset={offset1}). "
                                        f"There is no Substrand there.")
        if ss2 is None:
            raise IllegalDNADesignError(f"Cannot add half crossover at (helix={helix2}, offset={offset2}). "
                                        f"There is no Substrand there.")
        strand1 = ss1.strand()
        strand2 = ss2.strand()

        if strand1 == strand2:
            raise IllegalDNADesignError(f"Cannot add crossover from "
                                        f"(helix={helix1}, offset={offset1}) to "
                                        f"(helix={helix2}, offset={offset2}) "
                                        f"because that would join two Substrands "
                                        f"already on the same Strand! "
                                        f"Currently circular Strands are not supported. "
                                        f"Instead, try adding nicks first, or rearrange the order of "
                                        f"crossover addition, to ensure that all strands are "
                                        f"non-circular, even in intermediate stages.")

        if ss1.offset_3p() == offset1 and ss2.offset_5p() == offset2:
            strand_first = strand1
            strand_last = strand2
        elif ss1.offset_5p() == offset1 and ss2.offset_3p() == offset2:
            strand_first = strand2
            strand_last = strand1
        else:
            raise IllegalDNADesignError("Cannot add half crossover. Must have one substrand have its "
                                        "5' end at the given offset and the other with its 3' end at the "
                                        "given offset, but this is not the case.")

        new_substrands = strand_first.substrands + strand_last.substrands
        if strand_first.dna_sequence is None and strand_last.dna_sequence is None:
            new_dna = None
        elif strand_first.dna_sequence is not None and strand_last.dna_sequence is not None:
            new_dna = strand_first.dna_sequence + strand_last.dna_sequence
        else:
            raise IllegalDNADesignError('cannot add crossover between two strands if one has a DNA sequence '
                                        'and the other does not')
        new_strand = Strand(substrands=new_substrands, color=strand_first.color, dna_sequence=new_dna,
                            idt=strand_first.idt)

        self.strands.remove(strand_first)
        self.strands.remove(strand_last)
        self.strands.append(new_strand)

    def add_full_crossover(self, helix1: int, helix2: int, offset1: int, forward1: bool,
                           offset2: int = None, forward2: bool = None):
        """
        Adds two half-crossovers, one at `offset1` and another at `offset1`-1.
        Other arguments have the same meaning as in :py:meth:`DNADesign.add_half_crossover`.
        A nick is automatically added on helix `helix1` between
        `offset1` and `offset1`-1 if one is not already present,
        and similarly for `offset2` on helix `helix2`.
        """
        if offset2 is None:
            offset2 = offset1
        if forward2 is None:
            forward2 = not forward1
        for helix, forward, offset in [(helix1, forward1, offset1), (helix2, forward2, offset2)]:
            self._prepare_nicks_for_full_crossover(helix, forward, offset)
        self.add_half_crossover(helix1=helix1, helix2=helix2, offset1=offset1 - 1, offset2=offset2 - 1,
                                forward1=forward1, forward2=forward2)
        self.add_half_crossover(helix1=helix1, helix2=helix2, offset1=offset1, offset2=offset2,
                                forward1=forward1, forward2=forward2)

    def add_crossovers(self, crossovers: List[Crossover]):
        """
        Adds a list of :any:`Crossover`'s in batch.

        This helps to avoid problems where adding them one at a time
        creates an intermediate design with circular strands.

        :param crossovers: list of :any:`Crossover`'s to add. Its fields have the same meaning as in
            :py:meth:`DNADesign.add_half_crossover`
            and
            :py:meth:`DNADesign.add_full_crossover`,
            with the extra field `Crossover.half` indicating whether it represents a half or full crossover.

        """
        for crossover in crossovers:
            if not crossover.half:
                for helix, forward, offset in [(crossover.helix1, crossover.forward1, crossover.offset1),
                                               (crossover.helix2, crossover.forward2, crossover.offset2)]:
                    self._prepare_nicks_for_full_crossover(helix, forward, offset)

        for crossover in crossovers:
            if crossover.half:
                self.add_half_crossover(helix1=crossover.helix1, helix2=crossover.helix2,
                                        forward1=crossover.forward1, forward2=crossover.forward2,
                                        offset1=crossover.offset1, offset2=crossover.offset2)
            else:
                self.add_full_crossover(helix1=crossover.helix1, helix2=crossover.helix2,
                                        forward1=crossover.forward1, forward2=crossover.forward2,
                                        offset1=crossover.offset1, offset2=crossover.offset2)

    def _prepare_nicks_for_full_crossover(self, helix, forward, offset):
        substrand_right = self.substrand_at(helix, offset, forward)
        if substrand_right is None:
            raise IllegalDNADesignError(f'You tried to create a full crossover at '
                                        f'(helix={helix}, offset={offset}) '
                                        f'but there is no Strand there.')
        substrand_left = self.substrand_at(helix, offset - 1, forward)
        if substrand_left is None:
            raise IllegalDNADesignError(f'You tried to create a full crossover at '
                                        f'(helix={helix}, offset={offset}) '
                                        f'but there is no Strand at offset {offset - 1}.')
        if substrand_left == substrand_right:
            self.add_nick(helix, offset, forward)
        else:
            assert substrand_left.end == substrand_right.start

    def inline_deletions_insertions(self):
        """
        Converts deletions and insertions by "inlining" them. Insertions and deletions are removed,
        and their substrands have their lengths altered. Also, major tick marks on the helices will be
        shifted to preserve their adjacency to bases already present. For example, if there are major
        tick marks at 0, 8, 18, 24, and a deletion between 0 and 8, then
        the substrand is shortened by 1,
        the tick marks become 0, 7, 15, 23,
        and the helix's maximum offset is shrunk by 1.

        We assume that a major tick mark appears just to the LEFT of the offset it encodes,
        so the minimum and maximum offsets for tick marks are respectively the helix's minimum offset
        and 1 plus its maximum offset, the latter being just to the right of the last offset on the helix.
        """
        for helix in self.helices:
            self._inline_deletions_insertions_on_helix(helix)

    def _inline_deletions_insertions_on_helix(self, helix):
        ###################################################
        # first gather information before changing anything

        # gather all mods on helix
        deletions = [deletion for substrand in helix._substrands for deletion in substrand.deletions]
        insertions = [insertion for substrand in helix._substrands for insertion in substrand.insertions]

        # change max offset
        delta_length = sum(length for (offset, length) in insertions) - len(deletions)

        # combined collection of deletions/insertions into one dict mapping offset --> None/len, where
        # value of -1 indicates deletion, and otherwise is length of insertion
        dels_ins = dict()
        for deletion in deletions:
            dels_ins[deletion] = -1
        for insertion in insertions:
            dels_ins[insertion[0]] = insertion[1]

        # put offsets in sorted order
        dels_ins_offsets_sorted = sorted(dels_ins.keys())

        # fix helix major ticks
        major_ticks = sorted(helix.calculate_major_ticks(self.major_tick_distance))

        ###################################################
        # now that info is gathered, start changing things

        helix.max_offset += delta_length
        if len(major_ticks) > 0:
            major_tick_idx = 0
            delta_acc = 0  # accumulated delta; insertions add to this and deletions subtract from it
            for offset in dels_ins_offsets_sorted:
                # go to first major tick great than offset, updating passed ones by delta_acc
                while major_tick_idx < len(major_ticks) and major_ticks[major_tick_idx] <= offset:
                    major_ticks[major_tick_idx] += delta_acc
                    major_tick_idx += 1
                delta_acc += dels_ins[offset]
            # if necessary, update major ticks beyond last ins/del
            while major_tick_idx < len(major_ticks):
                major_ticks[major_tick_idx] += delta_acc
                major_tick_idx += 1
            # TODO: check if regularly spaced and reaching both ends, and if so set helix.major_tick_distance
            helix.major_ticks = major_ticks

        # fix substrand start/end offsets
        substrands = sorted(helix._substrands, key=lambda substrand: substrand.start)
        delta_acc = 0
        for substrand in substrands:
            substrand.start += delta_acc
            delta_acc += substrand.dna_length() - substrand.visual_length()
            substrand.end += delta_acc
            substrand.deletions = []
            substrand.insertions = []

    def reverse_all(self):
        """
        Reverses "polarity" of every :any:`Strand` in this :any:`DNADesign`.

        No attempt is made to make any assigned DNA sequences match by reversing or rearranging them.
        Every :any:`Strand` keeps the same DNA sequence it had before (unreversed), if one was assigned.
        It is recommended to assign/reassign DNA sequences *after* doing this operation.
        """
        for strand in self.strands:
            strand.reverse()


def _name_of_this_script() -> str:
    """Return name of the currently running script, WITHOUT the .py extension."""
    return os.path.basename(sys.argv[0])[:-3]


def _write_file_same_name_as_running_python_script(contents: str, extension: str, directory: str = '.',
                                                   filename=None):
    relative_filename = _get_filename_same_name_as_running_python_script(directory, extension, filename)
    with open(relative_filename, 'w') as out_file:
        out_file.write(contents)


def _get_filename_same_name_as_running_python_script(directory, extension, filename):
    if filename is None:
        filename = _name_of_this_script() + f'.{extension}'
    relative_filename = _create_directory_and_set_filename(directory, filename)
    return relative_filename


def _create_directory_and_set_filename(directory, filename):
    if not os.path.exists(directory):
        os.makedirs(directory)
    relative_filename = os.path.join(directory, filename)
    return relative_filename


def _remove_whitespace_and_uppercase(sequence):
    sequence = re.sub(r'\s*', '', sequence)
    sequence = sequence.upper()
    return sequence


def _pad_and_remove_whitespace(sequence, strand):
    sequence = _remove_whitespace_and_uppercase(sequence)
    padded_sequence = _pad_dna(sequence, strand.dna_length())
    return padded_sequence


def _pad_dna(sequence: str, length: int) -> str:
    """Return `sequence` modified to have length `length`.

    If len(sequence) < length, pad with  :py:data:`DNA_base_wildcard`.
    If len(sequence) > length, remove extra symbols."""
    if len(sequence) > length:
        sequence = sequence[:length]
    elif len(sequence) < length:
        sequence += DNA_base_wildcard * (length - len(sequence))
    return sequence


@dataclass
class Crossover:
    """
    A :any:`Crossover` object represents the parameters to the methods
    :py:meth:`DNADesign.add_half_crossover`
    and
    :py:meth:`DNADesign.add_full_crossover`,
    with one more field :py:data:`Crossover.half` to identify whether it is a half or full crossover.

    It is used in conjection with :py:meth:`DNADesign.add_crossovers` to add many crossovers in batch.
    This helps avoid the issue that adding crossovers one at a time can lead to an intermediate
    :any:`DNADesign` with circular strands, which are currently unsupported.
    """

    helix1: int
    """index of one helix of half crossover"""

    helix2: int
    """index of other helix of half crossover"""

    offset1: int
    """offset on `helix1` at which to add half crossover"""

    forward1: bool
    """direction of :any:`Strand` on `helix1` to which to add half crossover"""

    offset2: int = None
    """
    offset on `helix2` at which to add half crossover. If not specified, defaults to `offset1`
    """

    forward2: bool = None
    """
    direction of :any:`Strand` on `helix2` to which to add half crossover. 
    If not specified, defaults to the negation of `forward1`
    """

    half: bool = False
    """
    Indicates whether this is a half or full crossover.
    If not specified, defaults to ``False``.
    """

    def __post_init__(self):
        if self.offset2 is None:
            self.offset2 = self.offset1
        if self.forward2 is None:
            self.forward2 = not self.forward1


@dataclass
class DNAOrigamiDesign(DNADesign):
    """Subclass of :any:`DNADesign` that also defines a special "scaffold" strand as a field
    :py:data:`DNAOrigamiDesign.scaffold`.

    The field :py:data:`DNAOrigamiDesign.scaffold` can be set in the constructor,
    or it may be set after the :any:`DNAOrigamiDesign` is created.
    It should be assigned using the method :py:meth:`DNAOrigamiDesign.set_scaffold`.

    The :py:data:`Color` of the scaffold will be automatically assigned.
    To change from the default :any:`Color`,
    change the field :py:data:`Strand.color` of
    :py:data:`DNAOrigamiDesign.scaffold`
    *after* the :any:`DNAOrigamiDesign` is created.
    """

    scaffold: Strand = None
    """The scaffold :any:`Strand` of this :any:`DNAOrigamiDesign`.
    
    Must be an element of :py:data:`DNAOrigamiDesign.strands`."""

    def __post_init__(self):
        super().__post_init__()

        # XXX: it's not a great idea to allow scaffold to ne None, but this helps when someone wants
        # to create an origami by starting with several simple Strands and add nicks and crossovers.
        if self.scaffold is not None:
            self.scaffold.color = default_scaffold_color
            self.scaffold.is_scaffold = True
            if self.scaffold not in self.strands:
                raise StrandError(self.scaffold, 'scaffold strand not contained in DNAOrigamiDesigns.strands')

    def set_scaffold(self, scaffold: Strand):
        """
        Set the scaffold of this :any:`DNAOrigamiDesign`.

        :param scaffold: The scaffold :any:`Strand`.
        """
        self.scaffold = scaffold
        scaffold.is_scaffold = True
        self.scaffold.color = default_scaffold_color

    def assign_m13_to_scaffold(self):
        """
        Assigns the scaffold to be the sequence of M13: :py:data:`m13_sequence`.
        """
        self.assign_dna(self.scaffold, m13_sequence)

    # def to_json_serializable(self, suppress_indent=True):
    #     json_map = super().to_json_serializable(suppress_indent)
    #     json_map[is_origami_key] = True
    #     return json_map
