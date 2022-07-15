module_dict={}
module_dict["pyparsing/py.typed"]="""

"""
module_dict["pyparsing/unicode.py"]="""
# unicode.py

import sys
from itertools import filterfalse
from typing import List, Tuple, Union


class _lazyclassproperty:
    def __init__(self, fn):
        self.fn = fn
        self.__doc__ = fn.__doc__
        self.__name__ = fn.__name__

    def __get__(self, obj, cls):
        if cls is None:
            cls = type(obj)
        if not hasattr(cls, \"_intern\") or any(
            cls._intern is getattr(superclass, \"_intern\", [])
            for superclass in cls.__mro__[1:]
        ):
            cls._intern = {}
        attrname = self.fn.__name__
        if attrname not in cls._intern:
            cls._intern[attrname] = self.fn(cls)
        return cls._intern[attrname]


UnicodeRangeList = List[Union[Tuple[int, int], Tuple[int]]]


class unicode_set:
    \"\"\"
    A set of Unicode characters, for language-specific strings for
    ``alphas``, ``nums``, ``alphanums``, and ``printables``.
    A unicode_set is defined by a list of ranges in the Unicode character
    set, in a class attribute ``_ranges``. Ranges can be specified using
    2-tuples or a 1-tuple, such as::

        _ranges = [
            (0x0020, 0x007e),
            (0x00a0, 0x00ff),
            (0x0100,),
            ]

    Ranges are left- and right-inclusive. A 1-tuple of (x,) is treated as (x, x).

    A unicode set can also be defined using multiple inheritance of other unicode sets::

        class CJK(Chinese, Japanese, Korean):
            pass
    \"\"\"

    _ranges: UnicodeRangeList = []

    @_lazyclassproperty
    def _chars_for_ranges(cls):
        ret = []
        for cc in cls.__mro__:
            if cc is unicode_set:
                break
            for rr in getattr(cc, \"_ranges\", ()):
                ret.extend(range(rr[0], rr[-1] + 1))
        return [chr(c) for c in sorted(set(ret))]

    @_lazyclassproperty
    def printables(cls):
        \"all non-whitespace characters in this range\"
        return \"\".join(filterfalse(str.isspace, cls._chars_for_ranges))

    @_lazyclassproperty
    def alphas(cls):
        \"all alphabetic characters in this range\"
        return \"\".join(filter(str.isalpha, cls._chars_for_ranges))

    @_lazyclassproperty
    def nums(cls):
        \"all numeric digit characters in this range\"
        return \"\".join(filter(str.isdigit, cls._chars_for_ranges))

    @_lazyclassproperty
    def alphanums(cls):
        \"all alphanumeric characters in this range\"
        return cls.alphas + cls.nums

    @_lazyclassproperty
    def identchars(cls):
        \"all characters in this range that are valid identifier characters, plus underscore '_'\"
        return \"\".join(
            sorted(
                set(
                    \"\".join(filter(str.isidentifier, cls._chars_for_ranges))
                    + \"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzªµº\"
                    + \"ÀÁÂÃÄÅÆÇÈÉÊËÌÍÎÏÐÑÒÓÔÕÖØÙÚÛÜÝÞßàáâãäåæçèéêëìíîïðñòóôõöøùúûüýþÿ\"
                    + \"_\"
                )
            )
        )

    @_lazyclassproperty
    def identbodychars(cls):
        \"\"\"
        all characters in this range that are valid identifier body characters,
        plus the digits 0-9
        \"\"\"
        return \"\".join(
            sorted(
                set(
                    cls.identchars
                    + \"0123456789\"
                    + \"\".join(
                        [c for c in cls._chars_for_ranges if (\"_\" + c).isidentifier()]
                    )
                )
            )
        )


class pyparsing_unicode(unicode_set):
    \"\"\"
    A namespace class for defining common language unicode_sets.
    \"\"\"

    # fmt: off

    # define ranges in language character sets
    _ranges: UnicodeRangeList = [
        (0x0020, sys.maxunicode),
    ]

    class BasicMultilingualPlane(unicode_set):
        \"Unicode set for the Basic Multilingual Plane\"
        _ranges: UnicodeRangeList = [
            (0x0020, 0xFFFF),
        ]

    class Latin1(unicode_set):
        \"Unicode set for Latin-1 Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0020, 0x007E),
            (0x00A0, 0x00FF),
        ]

    class LatinA(unicode_set):
        \"Unicode set for Latin-A Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0100, 0x017F),
        ]

    class LatinB(unicode_set):
        \"Unicode set for Latin-B Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0180, 0x024F),
        ]

    class Greek(unicode_set):
        \"Unicode set for Greek Unicode Character Ranges\"
        _ranges: UnicodeRangeList = [
            (0x0342, 0x0345),
            (0x0370, 0x0377),
            (0x037A, 0x037F),
            (0x0384, 0x038A),
            (0x038C,),
            (0x038E, 0x03A1),
            (0x03A3, 0x03E1),
            (0x03F0, 0x03FF),
            (0x1D26, 0x1D2A),
            (0x1D5E,),
            (0x1D60,),
            (0x1D66, 0x1D6A),
            (0x1F00, 0x1F15),
            (0x1F18, 0x1F1D),
            (0x1F20, 0x1F45),
            (0x1F48, 0x1F4D),
            (0x1F50, 0x1F57),
            (0x1F59,),
            (0x1F5B,),
            (0x1F5D,),
            (0x1F5F, 0x1F7D),
            (0x1F80, 0x1FB4),
            (0x1FB6, 0x1FC4),
            (0x1FC6, 0x1FD3),
            (0x1FD6, 0x1FDB),
            (0x1FDD, 0x1FEF),
            (0x1FF2, 0x1FF4),
            (0x1FF6, 0x1FFE),
            (0x2129,),
            (0x2719, 0x271A),
            (0xAB65,),
            (0x10140, 0x1018D),
            (0x101A0,),
            (0x1D200, 0x1D245),
            (0x1F7A1, 0x1F7A7),
        ]

    class Cyrillic(unicode_set):
        \"Unicode set for Cyrillic Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0400, 0x052F),
            (0x1C80, 0x1C88),
            (0x1D2B,),
            (0x1D78,),
            (0x2DE0, 0x2DFF),
            (0xA640, 0xA672),
            (0xA674, 0xA69F),
            (0xFE2E, 0xFE2F),
        ]

    class Chinese(unicode_set):
        \"Unicode set for Chinese Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x2E80, 0x2E99),
            (0x2E9B, 0x2EF3),
            (0x31C0, 0x31E3),
            (0x3400, 0x4DB5),
            (0x4E00, 0x9FEF),
            (0xA700, 0xA707),
            (0xF900, 0xFA6D),
            (0xFA70, 0xFAD9),
            (0x16FE2, 0x16FE3),
            (0x1F210, 0x1F212),
            (0x1F214, 0x1F23B),
            (0x1F240, 0x1F248),
            (0x20000, 0x2A6D6),
            (0x2A700, 0x2B734),
            (0x2B740, 0x2B81D),
            (0x2B820, 0x2CEA1),
            (0x2CEB0, 0x2EBE0),
            (0x2F800, 0x2FA1D),
        ]

    class Japanese(unicode_set):
        \"Unicode set for Japanese Unicode Character Range, combining Kanji, Hiragana, and Katakana ranges\"
        _ranges: UnicodeRangeList = []

        class Kanji(unicode_set):
            \"Unicode set for Kanji Unicode Character Range\"
            _ranges: UnicodeRangeList = [
                (0x4E00, 0x9FBF),
                (0x3000, 0x303F),
            ]

        class Hiragana(unicode_set):
            \"Unicode set for Hiragana Unicode Character Range\"
            _ranges: UnicodeRangeList = [
                (0x3041, 0x3096),
                (0x3099, 0x30A0),
                (0x30FC,),
                (0xFF70,),
                (0x1B001,),
                (0x1B150, 0x1B152),
                (0x1F200,),
            ]

        class Katakana(unicode_set):
            \"Unicode set for Katakana  Unicode Character Range\"
            _ranges: UnicodeRangeList = [
                (0x3099, 0x309C),
                (0x30A0, 0x30FF),
                (0x31F0, 0x31FF),
                (0x32D0, 0x32FE),
                (0xFF65, 0xFF9F),
                (0x1B000,),
                (0x1B164, 0x1B167),
                (0x1F201, 0x1F202),
                (0x1F213,),
            ]

    class Hangul(unicode_set):
        \"Unicode set for Hangul (Korean) Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x1100, 0x11FF),
            (0x302E, 0x302F),
            (0x3131, 0x318E),
            (0x3200, 0x321C),
            (0x3260, 0x327B),
            (0x327E,),
            (0xA960, 0xA97C),
            (0xAC00, 0xD7A3),
            (0xD7B0, 0xD7C6),
            (0xD7CB, 0xD7FB),
            (0xFFA0, 0xFFBE),
            (0xFFC2, 0xFFC7),
            (0xFFCA, 0xFFCF),
            (0xFFD2, 0xFFD7),
            (0xFFDA, 0xFFDC),
        ]

    Korean = Hangul

    class CJK(Chinese, Japanese, Hangul):
        \"Unicode set for combined Chinese, Japanese, and Korean (CJK) Unicode Character Range\"

    class Thai(unicode_set):
        \"Unicode set for Thai Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0E01, 0x0E3A),
            (0x0E3F, 0x0E5B)
        ]

    class Arabic(unicode_set):
        \"Unicode set for Arabic Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0600, 0x061B),
            (0x061E, 0x06FF),
            (0x0700, 0x077F),
        ]

    class Hebrew(unicode_set):
        \"Unicode set for Hebrew Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0591, 0x05C7),
            (0x05D0, 0x05EA),
            (0x05EF, 0x05F4),
            (0xFB1D, 0xFB36),
            (0xFB38, 0xFB3C),
            (0xFB3E,),
            (0xFB40, 0xFB41),
            (0xFB43, 0xFB44),
            (0xFB46, 0xFB4F),
        ]

    class Devanagari(unicode_set):
        \"Unicode set for Devanagari Unicode Character Range\"
        _ranges: UnicodeRangeList = [
            (0x0900, 0x097F),
            (0xA8E0, 0xA8FF)
        ]

    # fmt: on


pyparsing_unicode.Japanese._ranges = (
    pyparsing_unicode.Japanese.Kanji._ranges
    + pyparsing_unicode.Japanese.Hiragana._ranges
    + pyparsing_unicode.Japanese.Katakana._ranges
)

pyparsing_unicode.BMP = pyparsing_unicode.BasicMultilingualPlane

# add language identifiers using language Unicode
pyparsing_unicode.العربية = pyparsing_unicode.Arabic
pyparsing_unicode.中文 = pyparsing_unicode.Chinese
pyparsing_unicode.кириллица = pyparsing_unicode.Cyrillic
pyparsing_unicode.Ελληνικά = pyparsing_unicode.Greek
pyparsing_unicode.עִברִית = pyparsing_unicode.Hebrew
pyparsing_unicode.日本語 = pyparsing_unicode.Japanese
pyparsing_unicode.Japanese.漢字 = pyparsing_unicode.Japanese.Kanji
pyparsing_unicode.Japanese.カタカナ = pyparsing_unicode.Japanese.Katakana
pyparsing_unicode.Japanese.ひらがな = pyparsing_unicode.Japanese.Hiragana
pyparsing_unicode.한국어 = pyparsing_unicode.Korean
pyparsing_unicode.ไทย = pyparsing_unicode.Thai
pyparsing_unicode.देवनागरी = pyparsing_unicode.Devanagari

"""
module_dict["pyparsing/core.py"]="""
#
# core.py
#
import os
import typing
from typing import (
    NamedTuple,
    Union,
    Callable,
    Any,
    Generator,
    Tuple,
    List,
    TextIO,
    Set,
    Sequence,
)
from abc import ABC, abstractmethod
from enum import Enum
import string
import copy
import warnings
import re
import sys
from collections.abc import Iterable
import traceback
import types
from operator import itemgetter
from functools import wraps
from threading import RLock
from pathlib import Path

from .util import (
    _FifoCache,
    _UnboundedCache,
    __config_flags,
    _collapse_string_to_ranges,
    _escape_regex_range_chars,
    _bslash,
    _flatten,
    LRUMemo as _LRUMemo,
    UnboundedMemo as _UnboundedMemo,
)
from .exceptions import *
from .actions import *
from .results import ParseResults, _ParseResultsWithOffset
from .unicode import pyparsing_unicode

_MAX_INT = sys.maxsize
str_type: Tuple[type, ...] = (str, bytes)

#
# Copyright (c) 2003-2022  Paul T. McGuire
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# \"Software\"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


if sys.version_info >= (3, 8):
    from functools import cached_property
else:

    class cached_property:
        def __init__(self, func):
            self._func = func

        def __get__(self, instance, owner=None):
            ret = instance.__dict__[self._func.__name__] = self._func(instance)
            return ret


class __compat__(__config_flags):
    \"\"\"
    A cross-version compatibility configuration for pyparsing features that will be
    released in a future version. By setting values in this configuration to True,
    those features can be enabled in prior versions for compatibility development
    and testing.

    - ``collect_all_And_tokens`` - flag to enable fix for Issue #63 that fixes erroneous grouping
      of results names when an :class:`And` expression is nested within an :class:`Or` or :class:`MatchFirst`;
      maintained for compatibility, but setting to ``False`` no longer restores pre-2.3.1
      behavior
    \"\"\"

    _type_desc = \"compatibility\"

    collect_all_And_tokens = True

    _all_names = [__ for __ in locals() if not __.startswith(\"_\")]
    _fixed_names = \"\"\"
        collect_all_And_tokens
        \"\"\".split()


class __diag__(__config_flags):
    _type_desc = \"diagnostic\"

    warn_multiple_tokens_in_named_alternation = False
    warn_ungrouped_named_tokens_in_collection = False
    warn_name_set_on_empty_Forward = False
    warn_on_parse_using_empty_Forward = False
    warn_on_assignment_to_Forward = False
    warn_on_multiple_string_args_to_oneof = False
    warn_on_match_first_with_lshift_operator = False
    enable_debug_on_named_expressions = False

    _all_names = [__ for __ in locals() if not __.startswith(\"_\")]
    _warning_names = [name for name in _all_names if name.startswith(\"warn\")]
    _debug_names = [name for name in _all_names if name.startswith(\"enable_debug\")]

    @classmethod
    def enable_all_warnings(cls) -> None:
        for name in cls._warning_names:
            cls.enable(name)


class Diagnostics(Enum):
    \"\"\"
    Diagnostic configuration (all default to disabled)
    - ``warn_multiple_tokens_in_named_alternation`` - flag to enable warnings when a results
      name is defined on a :class:`MatchFirst` or :class:`Or` expression with one or more :class:`And` subexpressions
    - ``warn_ungrouped_named_tokens_in_collection`` - flag to enable warnings when a results
      name is defined on a containing expression with ungrouped subexpressions that also
      have results names
    - ``warn_name_set_on_empty_Forward`` - flag to enable warnings when a :class:`Forward` is defined
      with a results name, but has no contents defined
    - ``warn_on_parse_using_empty_Forward`` - flag to enable warnings when a :class:`Forward` is
      defined in a grammar but has never had an expression attached to it
    - ``warn_on_assignment_to_Forward`` - flag to enable warnings when a :class:`Forward` is defined
      but is overwritten by assigning using ``'='`` instead of ``'<<='`` or ``'<<'``
    - ``warn_on_multiple_string_args_to_oneof`` - flag to enable warnings when :class:`one_of` is
      incorrectly called with multiple str arguments
    - ``enable_debug_on_named_expressions`` - flag to auto-enable debug on all subsequent
      calls to :class:`ParserElement.set_name`

    Diagnostics are enabled/disabled by calling :class:`enable_diag` and :class:`disable_diag`.
    All warnings can be enabled by calling :class:`enable_all_warnings`.
    \"\"\"

    warn_multiple_tokens_in_named_alternation = 0
    warn_ungrouped_named_tokens_in_collection = 1
    warn_name_set_on_empty_Forward = 2
    warn_on_parse_using_empty_Forward = 3
    warn_on_assignment_to_Forward = 4
    warn_on_multiple_string_args_to_oneof = 5
    warn_on_match_first_with_lshift_operator = 6
    enable_debug_on_named_expressions = 7


def enable_diag(diag_enum: Diagnostics) -> None:
    \"\"\"
    Enable a global pyparsing diagnostic flag (see :class:`Diagnostics`).
    \"\"\"
    __diag__.enable(diag_enum.name)


def disable_diag(diag_enum: Diagnostics) -> None:
    \"\"\"
    Disable a global pyparsing diagnostic flag (see :class:`Diagnostics`).
    \"\"\"
    __diag__.disable(diag_enum.name)


def enable_all_warnings() -> None:
    \"\"\"
    Enable all global pyparsing diagnostic warnings (see :class:`Diagnostics`).
    \"\"\"
    __diag__.enable_all_warnings()


# hide abstract class
del __config_flags


def _should_enable_warnings(
    cmd_line_warn_options: typing.Iterable[str], warn_env_var: typing.Optional[str]
) -> bool:
    enable = bool(warn_env_var)
    for warn_opt in cmd_line_warn_options:
        w_action, w_message, w_category, w_module, w_line = (warn_opt + \"::::\").split(
            \":\"
        )[:5]
        if not w_action.lower().startswith(\"i\") and (
            not (w_message or w_category or w_module) or w_module == \"pyparsing\"
        ):
            enable = True
        elif w_action.lower().startswith(\"i\") and w_module in (\"pyparsing\", \"\"):
            enable = False
    return enable


if _should_enable_warnings(
    sys.warnoptions, os.environ.get(\"PYPARSINGENABLEALLWARNINGS\")
):
    enable_all_warnings()


# build list of single arg builtins, that can be used as parse actions
_single_arg_builtins = {
    sum,
    len,
    sorted,
    reversed,
    list,
    tuple,
    set,
    any,
    all,
    min,
    max,
}

_generatorType = types.GeneratorType
ParseAction = Union[
    Callable[[], Any],
    Callable[[ParseResults], Any],
    Callable[[int, ParseResults], Any],
    Callable[[str, int, ParseResults], Any],
]
ParseCondition = Union[
    Callable[[], bool],
    Callable[[ParseResults], bool],
    Callable[[int, ParseResults], bool],
    Callable[[str, int, ParseResults], bool],
]
ParseFailAction = Callable[[str, int, \"ParserElement\", Exception], None]
DebugStartAction = Callable[[str, int, \"ParserElement\", bool], None]
DebugSuccessAction = Callable[
    [str, int, int, \"ParserElement\", ParseResults, bool], None
]
DebugExceptionAction = Callable[[str, int, \"ParserElement\", Exception, bool], None]


alphas = string.ascii_uppercase + string.ascii_lowercase
identchars = pyparsing_unicode.Latin1.identchars
identbodychars = pyparsing_unicode.Latin1.identbodychars
nums = \"0123456789\"
hexnums = nums + \"ABCDEFabcdef\"
alphanums = alphas + nums
printables = \"\".join([c for c in string.printable if c not in string.whitespace])

_trim_arity_call_line: traceback.StackSummary = None


def _trim_arity(func, max_limit=3):
    \"\"\"decorator to trim function calls to match the arity of the target\"\"\"
    global _trim_arity_call_line

    if func in _single_arg_builtins:
        return lambda s, l, t: func(t)

    limit = 0
    found_arity = False

    def extract_tb(tb, limit=0):
        frames = traceback.extract_tb(tb, limit=limit)
        frame_summary = frames[-1]
        return [frame_summary[:2]]

    # synthesize what would be returned by traceback.extract_stack at the call to
    # user's parse action 'func', so that we don't incur call penalty at parse time

    # fmt: off
    LINE_DIFF = 7
    # IF ANY CODE CHANGES, EVEN JUST COMMENTS OR BLANK LINES, BETWEEN THE NEXT LINE AND
    # THE CALL TO FUNC INSIDE WRAPPER, LINE_DIFF MUST BE MODIFIED!!!!
    _trim_arity_call_line = (_trim_arity_call_line or traceback.extract_stack(limit=2)[-1])
    pa_call_line_synth = (_trim_arity_call_line[0], _trim_arity_call_line[1] + LINE_DIFF)

    def wrapper(*args):
        nonlocal found_arity, limit
        while 1:
            try:
                ret = func(*args[limit:])
                found_arity = True
                return ret
            except TypeError as te:
                # re-raise TypeErrors if they did not come from our arity testing
                if found_arity:
                    raise
                else:
                    tb = te.__traceback__
                    trim_arity_type_error = (
                        extract_tb(tb, limit=2)[-1][:2] == pa_call_line_synth
                    )
                    del tb

                    if trim_arity_type_error:
                        if limit < max_limit:
                            limit += 1
                            continue

                    raise
    # fmt: on

    # copy func name to wrapper for sensible debug output
    # (can't use functools.wraps, since that messes with function signature)
    func_name = getattr(func, \"__name__\", getattr(func, \"__class__\").__name__)
    wrapper.__name__ = func_name
    wrapper.__doc__ = func.__doc__

    return wrapper


def condition_as_parse_action(
    fn: ParseCondition, message: str = None, fatal: bool = False
) -> ParseAction:
    \"\"\"
    Function to convert a simple predicate function that returns ``True`` or ``False``
    into a parse action. Can be used in places when a parse action is required
    and :class:`ParserElement.add_condition` cannot be used (such as when adding a condition
    to an operator level in :class:`infix_notation`).

    Optional keyword arguments:

    - ``message`` - define a custom message to be used in the raised exception
    - ``fatal`` - if True, will raise :class:`ParseFatalException` to stop parsing immediately;
      otherwise will raise :class:`ParseException`

    \"\"\"
    msg = message if message is not None else \"failed user-defined condition\"
    exc_type = ParseFatalException if fatal else ParseException
    fn = _trim_arity(fn)

    @wraps(fn)
    def pa(s, l, t):
        if not bool(fn(s, l, t)):
            raise exc_type(s, l, msg)

    return pa


def _default_start_debug_action(
    instring: str, loc: int, expr: \"ParserElement\", cache_hit: bool = False
):
    cache_hit_str = \"*\" if cache_hit else \"\"
    print(
        (
            \"{}Match {} at loc {}({},{})\\n  {}\\n  {}^\".format(
                cache_hit_str,
                expr,
                loc,
                lineno(loc, instring),
                col(loc, instring),
                line(loc, instring),
                \" \" * (col(loc, instring) - 1),
            )
        )
    )


def _default_success_debug_action(
    instring: str,
    startloc: int,
    endloc: int,
    expr: \"ParserElement\",
    toks: ParseResults,
    cache_hit: bool = False,
):
    cache_hit_str = \"*\" if cache_hit else \"\"
    print(\"{}Matched {} -> {}\".format(cache_hit_str, expr, toks.as_list()))


def _default_exception_debug_action(
    instring: str,
    loc: int,
    expr: \"ParserElement\",
    exc: Exception,
    cache_hit: bool = False,
):
    cache_hit_str = \"*\" if cache_hit else \"\"
    print(
        \"{}Match {} failed, {} raised: {}\".format(
            cache_hit_str, expr, type(exc).__name__, exc
        )
    )


def null_debug_action(*args):
    \"\"\"'Do-nothing' debug action, to suppress debugging output during parsing.\"\"\"


class ParserElement(ABC):
    \"\"\"Abstract base level parser element class.\"\"\"

    DEFAULT_WHITE_CHARS: str = \" \\n\\t\\r\"
    verbose_stacktrace: bool = False
    _literalStringClass: typing.Optional[type] = None

    @staticmethod
    def set_default_whitespace_chars(chars: str) -> None:
        r\"\"\"
        Overrides the default whitespace chars

        Example::

            # default whitespace chars are space, <TAB> and newline
            Word(alphas)[1, ...].parse_string(\"abc def\\nghi jkl\")  # -> ['abc', 'def', 'ghi', 'jkl']

            # change to just treat newline as significant
            ParserElement.set_default_whitespace_chars(\" \\t\")
            Word(alphas)[1, ...].parse_string(\"abc def\\nghi jkl\")  # -> ['abc', 'def']
        \"\"\"
        ParserElement.DEFAULT_WHITE_CHARS = chars

        # update whitespace all parse expressions defined in this module
        for expr in _builtin_exprs:
            if expr.copyDefaultWhiteChars:
                expr.whiteChars = set(chars)

    @staticmethod
    def inline_literals_using(cls: type) -> None:
        \"\"\"
        Set class to be used for inclusion of string literals into a parser.

        Example::

            # default literal class used is Literal
            integer = Word(nums)
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")

            date_str.parse_string(\"1999/12/31\")  # -> ['1999', '/', '12', '/', '31']


            # change to Suppress
            ParserElement.inline_literals_using(Suppress)
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")

            date_str.parse_string(\"1999/12/31\")  # -> ['1999', '12', '31']
        \"\"\"
        ParserElement._literalStringClass = cls

    class DebugActions(NamedTuple):
        debug_try: typing.Optional[DebugStartAction]
        debug_match: typing.Optional[DebugSuccessAction]
        debug_fail: typing.Optional[DebugExceptionAction]

    def __init__(self, savelist: bool = False):
        self.parseAction: List[ParseAction] = list()
        self.failAction: typing.Optional[ParseFailAction] = None
        self.customName = None
        self._defaultName = None
        self.resultsName = None
        self.saveAsList = savelist
        self.skipWhitespace = True
        self.whiteChars = set(ParserElement.DEFAULT_WHITE_CHARS)
        self.copyDefaultWhiteChars = True
        # used when checking for left-recursion
        self.mayReturnEmpty = False
        self.keepTabs = False
        self.ignoreExprs: List[\"ParserElement\"] = list()
        self.debug = False
        self.streamlined = False
        # optimize exception handling for subclasses that don't advance parse index
        self.mayIndexError = True
        self.errmsg = \"\"
        # mark results names as modal (report only last) or cumulative (list all)
        self.modalResults = True
        # custom debug actions
        self.debugActions = self.DebugActions(None, None, None)
        # avoid redundant calls to preParse
        self.callPreparse = True
        self.callDuringTry = False
        self.suppress_warnings_: List[Diagnostics] = []

    def suppress_warning(self, warning_type: Diagnostics) -> \"ParserElement\":
        \"\"\"
        Suppress warnings emitted for a particular diagnostic on this expression.

        Example::

            base = pp.Forward()
            base.suppress_warning(Diagnostics.warn_on_parse_using_empty_Forward)

            # statement would normally raise a warning, but is now suppressed
            print(base.parseString(\"x\"))

        \"\"\"
        self.suppress_warnings_.append(warning_type)
        return self

    def copy(self) -> \"ParserElement\":
        \"\"\"
        Make a copy of this :class:`ParserElement`.  Useful for defining
        different parse actions for the same parsing pattern, using copies of
        the original parse element.

        Example::

            integer = Word(nums).set_parse_action(lambda toks: int(toks[0]))
            integerK = integer.copy().add_parse_action(lambda toks: toks[0] * 1024) + Suppress(\"K\")
            integerM = integer.copy().add_parse_action(lambda toks: toks[0] * 1024 * 1024) + Suppress(\"M\")

            print((integerK | integerM | integer)[1, ...].parse_string(\"5K 100 640K 256M\"))

        prints::

            [5120, 100, 655360, 268435456]

        Equivalent form of ``expr.copy()`` is just ``expr()``::

            integerM = integer().add_parse_action(lambda toks: toks[0] * 1024 * 1024) + Suppress(\"M\")
        \"\"\"
        cpy = copy.copy(self)
        cpy.parseAction = self.parseAction[:]
        cpy.ignoreExprs = self.ignoreExprs[:]
        if self.copyDefaultWhiteChars:
            cpy.whiteChars = set(ParserElement.DEFAULT_WHITE_CHARS)
        return cpy

    def set_results_name(
        self, name: str, list_all_matches: bool = False, *, listAllMatches: bool = False
    ) -> \"ParserElement\":
        \"\"\"
        Define name for referencing matching tokens as a nested attribute
        of the returned parse results.

        Normally, results names are assigned as you would assign keys in a dict:
        any existing value is overwritten by later values. If it is necessary to
        keep all values captured for a particular results name, call ``set_results_name``
        with ``list_all_matches`` = True.

        NOTE: ``set_results_name`` returns a *copy* of the original :class:`ParserElement` object;
        this is so that the client can define a basic element, such as an
        integer, and reference it in multiple places with different names.

        You can also set results names using the abbreviated syntax,
        ``expr(\"name\")`` in place of ``expr.set_results_name(\"name\")``
        - see :class:`__call__`. If ``list_all_matches`` is required, use
        ``expr(\"name*\")``.

        Example::

            date_str = (integer.set_results_name(\"year\") + '/'
                        + integer.set_results_name(\"month\") + '/'
                        + integer.set_results_name(\"day\"))

            # equivalent form:
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")
        \"\"\"
        listAllMatches = listAllMatches or list_all_matches
        return self._setResultsName(name, listAllMatches)

    def _setResultsName(self, name, listAllMatches=False):
        if name is None:
            return self
        newself = self.copy()
        if name.endswith(\"*\"):
            name = name[:-1]
            listAllMatches = True
        newself.resultsName = name
        newself.modalResults = not listAllMatches
        return newself

    def set_break(self, break_flag: bool = True) -> \"ParserElement\":
        \"\"\"
        Method to invoke the Python pdb debugger when this element is
        about to be parsed. Set ``break_flag`` to ``True`` to enable, ``False`` to
        disable.
        \"\"\"
        if break_flag:
            _parseMethod = self._parse

            def breaker(instring, loc, doActions=True, callPreParse=True):
                import pdb

                # this call to pdb.set_trace() is intentional, not a checkin error
                pdb.set_trace()
                return _parseMethod(instring, loc, doActions, callPreParse)

            breaker._originalParseMethod = _parseMethod
            self._parse = breaker
        else:
            if hasattr(self._parse, \"_originalParseMethod\"):
                self._parse = self._parse._originalParseMethod
        return self

    def set_parse_action(self, *fns: ParseAction, **kwargs) -> \"ParserElement\":
        \"\"\"
        Define one or more actions to perform when successfully matching parse element definition.

        Parse actions can be called to perform data conversions, do extra validation,
        update external data structures, or enhance or replace the parsed tokens.
        Each parse action ``fn`` is a callable method with 0-3 arguments, called as
        ``fn(s, loc, toks)`` , ``fn(loc, toks)`` , ``fn(toks)`` , or just ``fn()`` , where:

        - s   = the original string being parsed (see note below)
        - loc = the location of the matching substring
        - toks = a list of the matched tokens, packaged as a :class:`ParseResults` object

        The parsed tokens are passed to the parse action as ParseResults. They can be
        modified in place using list-style append, extend, and pop operations to update
        the parsed list elements; and with dictionary-style item set and del operations
        to add, update, or remove any named results. If the tokens are modified in place,
        it is not necessary to return them with a return statement.

        Parse actions can also completely replace the given tokens, with another ``ParseResults``
        object, or with some entirely different object (common for parse actions that perform data
        conversions). A convenient way to build a new parse result is to define the values
        using a dict, and then create the return value using :class:`ParseResults.from_dict`.

        If None is passed as the ``fn`` parse action, all previously added parse actions for this
        expression are cleared.

        Optional keyword arguments:

        - call_during_try = (default= ``False``) indicate if parse action should be run during
          lookaheads and alternate testing. For parse actions that have side effects, it is
          important to only call the parse action once it is determined that it is being
          called as part of a successful parse. For parse actions that perform additional
          validation, then call_during_try should be passed as True, so that the validation
          code is included in the preliminary \"try\" parses.

        Note: the default parsing behavior is to expand tabs in the input string
        before starting the parsing process.  See :class:`parse_string` for more
        information on parsing strings containing ``<TAB>`` s, and suggested
        methods to maintain a consistent view of the parsed string, the parse
        location, and line and column positions within the parsed string.

        Example::

            # parse dates in the form YYYY/MM/DD

            # use parse action to convert toks from str to int at parse time
            def convert_to_int(toks):
                return int(toks[0])

            # use a parse action to verify that the date is a valid date
            def is_valid_date(instring, loc, toks):
                from datetime import date
                year, month, day = toks[::2]
                try:
                    date(year, month, day)
                except ValueError:
                    raise ParseException(instring, loc, \"invalid date given\")

            integer = Word(nums)
            date_str = integer + '/' + integer + '/' + integer

            # add parse actions
            integer.set_parse_action(convert_to_int)
            date_str.set_parse_action(is_valid_date)

            # note that integer fields are now ints, not strings
            date_str.run_tests('''
                # successful parse - note that integer fields were converted to ints
                1999/12/31

                # fail - invalid date
                1999/13/31
                ''')
        \"\"\"
        if list(fns) == [None]:
            self.parseAction = []
        else:
            if not all(callable(fn) for fn in fns):
                raise TypeError(\"parse actions must be callable\")
            self.parseAction = [_trim_arity(fn) for fn in fns]
            self.callDuringTry = kwargs.get(
                \"call_during_try\", kwargs.get(\"callDuringTry\", False)
            )
        return self

    def add_parse_action(self, *fns: ParseAction, **kwargs) -> \"ParserElement\":
        \"\"\"
        Add one or more parse actions to expression's list of parse actions. See :class:`set_parse_action`.

        See examples in :class:`copy`.
        \"\"\"
        self.parseAction += [_trim_arity(fn) for fn in fns]
        self.callDuringTry = self.callDuringTry or kwargs.get(
            \"call_during_try\", kwargs.get(\"callDuringTry\", False)
        )
        return self

    def add_condition(self, *fns: ParseCondition, **kwargs) -> \"ParserElement\":
        \"\"\"Add a boolean predicate function to expression's list of parse actions. See
        :class:`set_parse_action` for function call signatures. Unlike ``set_parse_action``,
        functions passed to ``add_condition`` need to return boolean success/fail of the condition.

        Optional keyword arguments:

        - message = define a custom message to be used in the raised exception
        - fatal = if True, will raise ParseFatalException to stop parsing immediately; otherwise will raise
          ParseException
        - call_during_try = boolean to indicate if this method should be called during internal tryParse calls,
          default=False

        Example::

            integer = Word(nums).set_parse_action(lambda toks: int(toks[0]))
            year_int = integer.copy()
            year_int.add_condition(lambda toks: toks[0] >= 2000, message=\"Only support years 2000 and later\")
            date_str = year_int + '/' + integer + '/' + integer

            result = date_str.parse_string(\"1999/12/31\")  # -> Exception: Only support years 2000 and later (at char 0),
                                                                         (line:1, col:1)
        \"\"\"
        for fn in fns:
            self.parseAction.append(
                condition_as_parse_action(
                    fn, message=kwargs.get(\"message\"), fatal=kwargs.get(\"fatal\", False)
                )
            )

        self.callDuringTry = self.callDuringTry or kwargs.get(
            \"call_during_try\", kwargs.get(\"callDuringTry\", False)
        )
        return self

    def set_fail_action(self, fn: ParseFailAction) -> \"ParserElement\":
        \"\"\"
        Define action to perform if parsing fails at this expression.
        Fail acton fn is a callable function that takes the arguments
        ``fn(s, loc, expr, err)`` where:

        - s = string being parsed
        - loc = location where expression match was attempted and failed
        - expr = the parse expression that failed
        - err = the exception thrown

        The function returns no value.  It may throw :class:`ParseFatalException`
        if it is desired to stop parsing immediately.\"\"\"
        self.failAction = fn
        return self

    def _skipIgnorables(self, instring, loc):
        exprsFound = True
        while exprsFound:
            exprsFound = False
            for e in self.ignoreExprs:
                try:
                    while 1:
                        loc, dummy = e._parse(instring, loc)
                        exprsFound = True
                except ParseException:
                    pass
        return loc

    def preParse(self, instring, loc):
        if self.ignoreExprs:
            loc = self._skipIgnorables(instring, loc)

        if self.skipWhitespace:
            instrlen = len(instring)
            white_chars = self.whiteChars
            while loc < instrlen and instring[loc] in white_chars:
                loc += 1

        return loc

    def parseImpl(self, instring, loc, doActions=True):
        return loc, []

    def postParse(self, instring, loc, tokenlist):
        return tokenlist

    # @profile
    def _parseNoCache(
        self, instring, loc, doActions=True, callPreParse=True
    ) -> Tuple[int, ParseResults]:
        TRY, MATCH, FAIL = 0, 1, 2
        debugging = self.debug  # and doActions)
        len_instring = len(instring)

        if debugging or self.failAction:
            # print(\"Match {} at loc {}({}, {})\".format(self, loc, lineno(loc, instring), col(loc, instring)))
            try:
                if callPreParse and self.callPreparse:
                    pre_loc = self.preParse(instring, loc)
                else:
                    pre_loc = loc
                tokens_start = pre_loc
                if self.debugActions.debug_try:
                    self.debugActions.debug_try(instring, tokens_start, self, False)
                if self.mayIndexError or pre_loc >= len_instring:
                    try:
                        loc, tokens = self.parseImpl(instring, pre_loc, doActions)
                    except IndexError:
                        raise ParseException(instring, len_instring, self.errmsg, self)
                else:
                    loc, tokens = self.parseImpl(instring, pre_loc, doActions)
            except Exception as err:
                # print(\"Exception raised:\", err)
                if self.debugActions.debug_fail:
                    self.debugActions.debug_fail(
                        instring, tokens_start, self, err, False
                    )
                if self.failAction:
                    self.failAction(instring, tokens_start, self, err)
                raise
        else:
            if callPreParse and self.callPreparse:
                pre_loc = self.preParse(instring, loc)
            else:
                pre_loc = loc
            tokens_start = pre_loc
            if self.mayIndexError or pre_loc >= len_instring:
                try:
                    loc, tokens = self.parseImpl(instring, pre_loc, doActions)
                except IndexError:
                    raise ParseException(instring, len_instring, self.errmsg, self)
            else:
                loc, tokens = self.parseImpl(instring, pre_loc, doActions)

        tokens = self.postParse(instring, loc, tokens)

        ret_tokens = ParseResults(
            tokens, self.resultsName, asList=self.saveAsList, modal=self.modalResults
        )
        if self.parseAction and (doActions or self.callDuringTry):
            if debugging:
                try:
                    for fn in self.parseAction:
                        try:
                            tokens = fn(instring, tokens_start, ret_tokens)
                        except IndexError as parse_action_exc:
                            exc = ParseException(\"exception raised in parse action\")
                            raise exc from parse_action_exc

                        if tokens is not None and tokens is not ret_tokens:
                            ret_tokens = ParseResults(
                                tokens,
                                self.resultsName,
                                asList=self.saveAsList
                                and isinstance(tokens, (ParseResults, list)),
                                modal=self.modalResults,
                            )
                except Exception as err:
                    # print \"Exception raised in user parse action:\", err
                    if self.debugActions.debug_fail:
                        self.debugActions.debug_fail(
                            instring, tokens_start, self, err, False
                        )
                    raise
            else:
                for fn in self.parseAction:
                    try:
                        tokens = fn(instring, tokens_start, ret_tokens)
                    except IndexError as parse_action_exc:
                        exc = ParseException(\"exception raised in parse action\")
                        raise exc from parse_action_exc

                    if tokens is not None and tokens is not ret_tokens:
                        ret_tokens = ParseResults(
                            tokens,
                            self.resultsName,
                            asList=self.saveAsList
                            and isinstance(tokens, (ParseResults, list)),
                            modal=self.modalResults,
                        )
        if debugging:
            # print(\"Matched\", self, \"->\", ret_tokens.as_list())
            if self.debugActions.debug_match:
                self.debugActions.debug_match(
                    instring, tokens_start, loc, self, ret_tokens, False
                )

        return loc, ret_tokens

    def try_parse(self, instring: str, loc: int, raise_fatal: bool = False) -> int:
        try:
            return self._parse(instring, loc, doActions=False)[0]
        except ParseFatalException:
            if raise_fatal:
                raise
            raise ParseException(instring, loc, self.errmsg, self)

    def can_parse_next(self, instring: str, loc: int) -> bool:
        try:
            self.try_parse(instring, loc)
        except (ParseException, IndexError):
            return False
        else:
            return True

    # cache for left-recursion in Forward references
    recursion_lock = RLock()
    recursion_memos: typing.Dict[
        Tuple[int, \"Forward\", bool], Tuple[int, Union[ParseResults, Exception]]
    ] = {}

    # argument cache for optimizing repeated calls when backtracking through recursive expressions
    packrat_cache = (
        {}
    )  # this is set later by enabled_packrat(); this is here so that reset_cache() doesn't fail
    packrat_cache_lock = RLock()
    packrat_cache_stats = [0, 0]

    # this method gets repeatedly called during backtracking with the same arguments -
    # we can cache these arguments and save ourselves the trouble of re-parsing the contained expression
    def _parseCache(
        self, instring, loc, doActions=True, callPreParse=True
    ) -> Tuple[int, ParseResults]:
        HIT, MISS = 0, 1
        TRY, MATCH, FAIL = 0, 1, 2
        lookup = (self, instring, loc, callPreParse, doActions)
        with ParserElement.packrat_cache_lock:
            cache = ParserElement.packrat_cache
            value = cache.get(lookup)
            if value is cache.not_in_cache:
                ParserElement.packrat_cache_stats[MISS] += 1
                try:
                    value = self._parseNoCache(instring, loc, doActions, callPreParse)
                except ParseBaseException as pe:
                    # cache a copy of the exception, without the traceback
                    cache.set(lookup, pe.__class__(*pe.args))
                    raise
                else:
                    cache.set(lookup, (value[0], value[1].copy(), loc))
                    return value
            else:
                ParserElement.packrat_cache_stats[HIT] += 1
                if self.debug and self.debugActions.debug_try:
                    try:
                        self.debugActions.debug_try(instring, loc, self, cache_hit=True)
                    except TypeError:
                        pass
                if isinstance(value, Exception):
                    if self.debug and self.debugActions.debug_fail:
                        try:
                            self.debugActions.debug_fail(
                                instring, loc, self, value, cache_hit=True
                            )
                        except TypeError:
                            pass
                    raise value

                loc_, result, endloc = value[0], value[1].copy(), value[2]
                if self.debug and self.debugActions.debug_match:
                    try:
                        self.debugActions.debug_match(
                            instring, loc_, endloc, self, result, cache_hit=True
                        )
                    except TypeError:
                        pass

                return loc_, result

    _parse = _parseNoCache

    @staticmethod
    def reset_cache() -> None:
        ParserElement.packrat_cache.clear()
        ParserElement.packrat_cache_stats[:] = [0] * len(
            ParserElement.packrat_cache_stats
        )
        ParserElement.recursion_memos.clear()

    _packratEnabled = False
    _left_recursion_enabled = False

    @staticmethod
    def disable_memoization() -> None:
        \"\"\"
        Disables active Packrat or Left Recursion parsing and their memoization

        This method also works if neither Packrat nor Left Recursion are enabled.
        This makes it safe to call before activating Packrat nor Left Recursion
        to clear any previous settings.
        \"\"\"
        ParserElement.reset_cache()
        ParserElement._left_recursion_enabled = False
        ParserElement._packratEnabled = False
        ParserElement._parse = ParserElement._parseNoCache

    @staticmethod
    def enable_left_recursion(
        cache_size_limit: typing.Optional[int] = None, *, force=False
    ) -> None:
        \"\"\"
        Enables \"bounded recursion\" parsing, which allows for both direct and indirect
        left-recursion. During parsing, left-recursive :class:`Forward` elements are
        repeatedly matched with a fixed recursion depth that is gradually increased
        until finding the longest match.

        Example::

            import pyparsing as pp
            pp.ParserElement.enable_left_recursion()

            E = pp.Forward(\"E\")
            num = pp.Word(pp.nums)
            # match `num`, or `num '+' num`, or `num '+' num '+' num`, ...
            E <<= E + '+' - num | num

            print(E.parse_string(\"1+2+3\"))

        Recursion search naturally memoizes matches of ``Forward`` elements and may
        thus skip reevaluation of parse actions during backtracking. This may break
        programs with parse actions which rely on strict ordering of side-effects.

        Parameters:

        - cache_size_limit - (default=``None``) - memoize at most this many
          ``Forward`` elements during matching; if ``None`` (the default),
          memoize all ``Forward`` elements.

        Bounded Recursion parsing works similar but not identical to Packrat parsing,
        thus the two cannot be used together. Use ``force=True`` to disable any
        previous, conflicting settings.
        \"\"\"
        if force:
            ParserElement.disable_memoization()
        elif ParserElement._packratEnabled:
            raise RuntimeError(\"Packrat and Bounded Recursion are not compatible\")
        if cache_size_limit is None:
            ParserElement.recursion_memos = _UnboundedMemo()
        elif cache_size_limit > 0:
            ParserElement.recursion_memos = _LRUMemo(capacity=cache_size_limit)
        else:
            raise NotImplementedError(\"Memo size of %s\" % cache_size_limit)
        ParserElement._left_recursion_enabled = True

    @staticmethod
    def enable_packrat(cache_size_limit: int = 128, *, force: bool = False) -> None:
        \"\"\"
        Enables \"packrat\" parsing, which adds memoizing to the parsing logic.
        Repeated parse attempts at the same string location (which happens
        often in many complex grammars) can immediately return a cached value,
        instead of re-executing parsing/validating code.  Memoizing is done of
        both valid results and parsing exceptions.

        Parameters:

        - cache_size_limit - (default= ``128``) - if an integer value is provided
          will limit the size of the packrat cache; if None is passed, then
          the cache size will be unbounded; if 0 is passed, the cache will
          be effectively disabled.

        This speedup may break existing programs that use parse actions that
        have side-effects.  For this reason, packrat parsing is disabled when
        you first import pyparsing.  To activate the packrat feature, your
        program must call the class method :class:`ParserElement.enable_packrat`.
        For best results, call ``enable_packrat()`` immediately after
        importing pyparsing.

        Example::

            import pyparsing
            pyparsing.ParserElement.enable_packrat()

        Packrat parsing works similar but not identical to Bounded Recursion parsing,
        thus the two cannot be used together. Use ``force=True`` to disable any
        previous, conflicting settings.
        \"\"\"
        if force:
            ParserElement.disable_memoization()
        elif ParserElement._left_recursion_enabled:
            raise RuntimeError(\"Packrat and Bounded Recursion are not compatible\")
        if not ParserElement._packratEnabled:
            ParserElement._packratEnabled = True
            if cache_size_limit is None:
                ParserElement.packrat_cache = _UnboundedCache()
            else:
                ParserElement.packrat_cache = _FifoCache(cache_size_limit)
            ParserElement._parse = ParserElement._parseCache

    def parse_string(
        self, instring: str, parse_all: bool = False, *, parseAll: bool = False
    ) -> ParseResults:
        \"\"\"
        Parse a string with respect to the parser definition. This function is intended as the primary interface to the
        client code.

        :param instring: The input string to be parsed.
        :param parse_all: If set, the entire input string must match the grammar.
        :param parseAll: retained for pre-PEP8 compatibility, will be removed in a future release.
        :raises ParseException: Raised if ``parse_all`` is set and the input string does not match the whole grammar.
        :returns: the parsed data as a :class:`ParseResults` object, which may be accessed as a `list`, a `dict`, or
          an object with attributes if the given parser includes results names.

        If the input string is required to match the entire grammar, ``parse_all`` flag must be set to ``True``. This
        is also equivalent to ending the grammar with :class:`StringEnd`().

        To report proper column numbers, ``parse_string`` operates on a copy of the input string where all tabs are
        converted to spaces (8 spaces per tab, as per the default in ``string.expandtabs``). If the input string
        contains tabs and the grammar uses parse actions that use the ``loc`` argument to index into the string
        being parsed, one can ensure a consistent view of the input string by doing one of the following:

        - calling ``parse_with_tabs`` on your grammar before calling ``parse_string`` (see :class:`parse_with_tabs`),
        - define your parse action using the full ``(s,loc,toks)`` signature, and reference the input string using the
          parse action's ``s`` argument, or
        - explicitly expand the tabs in your input string before calling ``parse_string``.

        Examples:

        By default, partial matches are OK.

        >>> res = Word('a').parse_string('aaaaabaaa')
        >>> print(res)
        ['aaaaa']

        The parsing behavior varies by the inheriting class of this abstract class. Please refer to the children
        directly to see more examples.

        It raises an exception if parse_all flag is set and instring does not match the whole grammar.

        >>> res = Word('a').parse_string('aaaaabaaa', parse_all=True)
        Traceback (most recent call last):
        ...
        pyparsing.ParseException: Expected end of text, found 'b'  (at char 5), (line:1, col:6)
        \"\"\"
        parseAll = parse_all or parseAll

        ParserElement.reset_cache()
        if not self.streamlined:
            self.streamline()
        for e in self.ignoreExprs:
            e.streamline()
        if not self.keepTabs:
            instring = instring.expandtabs()
        try:
            loc, tokens = self._parse(instring, 0)
            if parseAll:
                loc = self.preParse(instring, loc)
                se = Empty() + StringEnd()
                se._parse(instring, loc)
        except ParseBaseException as exc:
            if ParserElement.verbose_stacktrace:
                raise
            else:
                # catch and re-raise exception from here, clearing out pyparsing internal stack trace
                raise exc.with_traceback(None)
        else:
            return tokens

    def scan_string(
        self,
        instring: str,
        max_matches: int = _MAX_INT,
        overlap: bool = False,
        *,
        debug: bool = False,
        maxMatches: int = _MAX_INT,
    ) -> Generator[Tuple[ParseResults, int, int], None, None]:
        \"\"\"
        Scan the input string for expression matches.  Each match will return the
        matching tokens, start location, and end location.  May be called with optional
        ``max_matches`` argument, to clip scanning after 'n' matches are found.  If
        ``overlap`` is specified, then overlapping matches will be reported.

        Note that the start and end locations are reported relative to the string
        being parsed.  See :class:`parse_string` for more information on parsing
        strings with embedded tabs.

        Example::

            source = \"sldjf123lsdjjkf345sldkjf879lkjsfd987\"
            print(source)
            for tokens, start, end in Word(alphas).scan_string(source):
                print(' '*start + '^'*(end-start))
                print(' '*start + tokens[0])

        prints::

            sldjf123lsdjjkf345sldkjf879lkjsfd987
            ^^^^^
            sldjf
                    ^^^^^^^
                    lsdjjkf
                              ^^^^^^
                              sldkjf
                                       ^^^^^^
                                       lkjsfd
        \"\"\"
        maxMatches = min(maxMatches, max_matches)
        if not self.streamlined:
            self.streamline()
        for e in self.ignoreExprs:
            e.streamline()

        if not self.keepTabs:
            instring = str(instring).expandtabs()
        instrlen = len(instring)
        loc = 0
        preparseFn = self.preParse
        parseFn = self._parse
        ParserElement.resetCache()
        matches = 0
        try:
            while loc <= instrlen and matches < maxMatches:
                try:
                    preloc = preparseFn(instring, loc)
                    nextLoc, tokens = parseFn(instring, preloc, callPreParse=False)
                except ParseException:
                    loc = preloc + 1
                else:
                    if nextLoc > loc:
                        matches += 1
                        if debug:
                            print(
                                {
                                    \"tokens\": tokens.asList(),
                                    \"start\": preloc,
                                    \"end\": nextLoc,
                                }
                            )
                        yield tokens, preloc, nextLoc
                        if overlap:
                            nextloc = preparseFn(instring, loc)
                            if nextloc > loc:
                                loc = nextLoc
                            else:
                                loc += 1
                        else:
                            loc = nextLoc
                    else:
                        loc = preloc + 1
        except ParseBaseException as exc:
            if ParserElement.verbose_stacktrace:
                raise
            else:
                # catch and re-raise exception from here, clears out pyparsing internal stack trace
                raise exc.with_traceback(None)

    def transform_string(self, instring: str, *, debug: bool = False) -> str:
        \"\"\"
        Extension to :class:`scan_string`, to modify matching text with modified tokens that may
        be returned from a parse action.  To use ``transform_string``, define a grammar and
        attach a parse action to it that modifies the returned token list.
        Invoking ``transform_string()`` on a target string will then scan for matches,
        and replace the matched text patterns according to the logic in the parse
        action.  ``transform_string()`` returns the resulting transformed string.

        Example::

            wd = Word(alphas)
            wd.set_parse_action(lambda toks: toks[0].title())

            print(wd.transform_string(\"now is the winter of our discontent made glorious summer by this sun of york.\"))

        prints::

            Now Is The Winter Of Our Discontent Made Glorious Summer By This Sun Of York.
        \"\"\"
        out: List[str] = []
        lastE = 0
        # force preservation of <TAB>s, to minimize unwanted transformation of string, and to
        # keep string locs straight between transform_string and scan_string
        self.keepTabs = True
        try:
            for t, s, e in self.scan_string(instring, debug=debug):
                out.append(instring[lastE:s])
                if t:
                    if isinstance(t, ParseResults):
                        out += t.as_list()
                    elif isinstance(t, Iterable) and not isinstance(t, str_type):
                        out.extend(t)
                    else:
                        out.append(t)
                lastE = e
            out.append(instring[lastE:])
            out = [o for o in out if o]
            return \"\".join([str(s) for s in _flatten(out)])
        except ParseBaseException as exc:
            if ParserElement.verbose_stacktrace:
                raise
            else:
                # catch and re-raise exception from here, clears out pyparsing internal stack trace
                raise exc.with_traceback(None)

    def search_string(
        self,
        instring: str,
        max_matches: int = _MAX_INT,
        *,
        debug: bool = False,
        maxMatches: int = _MAX_INT,
    ) -> ParseResults:
        \"\"\"
        Another extension to :class:`scan_string`, simplifying the access to the tokens found
        to match the given parse expression.  May be called with optional
        ``max_matches`` argument, to clip searching after 'n' matches are found.

        Example::

            # a capitalized word starts with an uppercase letter, followed by zero or more lowercase letters
            cap_word = Word(alphas.upper(), alphas.lower())

            print(cap_word.search_string(\"More than Iron, more than Lead, more than Gold I need Electricity\"))

            # the sum() builtin can be used to merge results into a single ParseResults object
            print(sum(cap_word.search_string(\"More than Iron, more than Lead, more than Gold I need Electricity\")))

        prints::

            [['More'], ['Iron'], ['Lead'], ['Gold'], ['I'], ['Electricity']]
            ['More', 'Iron', 'Lead', 'Gold', 'I', 'Electricity']
        \"\"\"
        maxMatches = min(maxMatches, max_matches)
        try:
            return ParseResults(
                [t for t, s, e in self.scan_string(instring, maxMatches, debug=debug)]
            )
        except ParseBaseException as exc:
            if ParserElement.verbose_stacktrace:
                raise
            else:
                # catch and re-raise exception from here, clears out pyparsing internal stack trace
                raise exc.with_traceback(None)

    def split(
        self,
        instring: str,
        maxsplit: int = _MAX_INT,
        include_separators: bool = False,
        *,
        includeSeparators=False,
    ) -> Generator[str, None, None]:
        \"\"\"
        Generator method to split a string using the given expression as a separator.
        May be called with optional ``maxsplit`` argument, to limit the number of splits;
        and the optional ``include_separators`` argument (default= ``False``), if the separating
        matching text should be included in the split results.

        Example::

            punc = one_of(list(\".,;:/-!?\"))
            print(list(punc.split(\"This, this?, this sentence, is badly punctuated!\")))

        prints::

            ['This', ' this', '', ' this sentence', ' is badly punctuated', '']
        \"\"\"
        includeSeparators = includeSeparators or include_separators
        last = 0
        for t, s, e in self.scan_string(instring, max_matches=maxsplit):
            yield instring[last:s]
            if includeSeparators:
                yield t[0]
            last = e
        yield instring[last:]

    def __add__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``+`` operator - returns :class:`And`. Adding strings to a :class:`ParserElement`
        converts them to :class:`Literal`s by default.

        Example::

            greet = Word(alphas) + \",\" + Word(alphas) + \"!\"
            hello = \"Hello, World!\"
            print(hello, \"->\", greet.parse_string(hello))

        prints::

            Hello, World! -> ['Hello', ',', 'World', '!']

        ``...`` may be used as a parse expression as a short form of :class:`SkipTo`.

            Literal('start') + ... + Literal('end')

        is equivalent to:

            Literal('start') + SkipTo('end')(\"_skipped*\") + Literal('end')

        Note that the skipped text is returned with '_skipped' as a results name,
        and to support having multiple skips in the same parser, the value returned is
        a list of all skipped text.
        \"\"\"
        if other is Ellipsis:
            return _PendingSkip(self)

        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return And([self, other])

    def __radd__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``+`` operator when left operand is not a :class:`ParserElement`
        \"\"\"
        if other is Ellipsis:
            return SkipTo(self)(\"_skipped*\") + self

        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return other + self

    def __sub__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``-`` operator, returns :class:`And` with error stop
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return self + And._ErrorStop() + other

    def __rsub__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``-`` operator when left operand is not a :class:`ParserElement`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return other - self

    def __mul__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``*`` operator, allows use of ``expr * 3`` in place of
        ``expr + expr + expr``.  Expressions may also be multiplied by a 2-integer
        tuple, similar to ``{min, max}`` multipliers in regular expressions.  Tuples
        may also include ``None`` as in:
        - ``expr*(n, None)`` or ``expr*(n, )`` is equivalent
             to ``expr*n + ZeroOrMore(expr)``
             (read as \"at least n instances of ``expr``\")
        - ``expr*(None, n)`` is equivalent to ``expr*(0, n)``
             (read as \"0 to n instances of ``expr``\")
        - ``expr*(None, None)`` is equivalent to ``ZeroOrMore(expr)``
        - ``expr*(1, None)`` is equivalent to ``OneOrMore(expr)``

        Note that ``expr*(None, n)`` does not raise an exception if
        more than n exprs exist in the input stream; that is,
        ``expr*(None, n)`` does not enforce a maximum number of expr
        occurrences.  If this behavior is desired, then write
        ``expr*(None, n) + ~expr``
        \"\"\"
        if other is Ellipsis:
            other = (0, None)
        elif isinstance(other, tuple) and other[:1] == (Ellipsis,):
            other = ((0,) + other[1:] + (None,))[:2]

        if isinstance(other, int):
            minElements, optElements = other, 0
        elif isinstance(other, tuple):
            other = tuple(o if o is not Ellipsis else None for o in other)
            other = (other + (None, None))[:2]
            if other[0] is None:
                other = (0, other[1])
            if isinstance(other[0], int) and other[1] is None:
                if other[0] == 0:
                    return ZeroOrMore(self)
                if other[0] == 1:
                    return OneOrMore(self)
                else:
                    return self * other[0] + ZeroOrMore(self)
            elif isinstance(other[0], int) and isinstance(other[1], int):
                minElements, optElements = other
                optElements -= minElements
            else:
                raise TypeError(
                    \"cannot multiply ParserElement and ({}) objects\".format(
                        \",\".join(type(item).__name__ for item in other)
                    )
                )
        else:
            raise TypeError(
                \"cannot multiply ParserElement and {} objects\".format(
                    type(other).__name__
                )
            )

        if minElements < 0:
            raise ValueError(\"cannot multiply ParserElement by negative value\")
        if optElements < 0:
            raise ValueError(
                \"second tuple value must be greater or equal to first tuple value\"
            )
        if minElements == optElements == 0:
            return And([])

        if optElements:

            def makeOptionalList(n):
                if n > 1:
                    return Opt(self + makeOptionalList(n - 1))
                else:
                    return Opt(self)

            if minElements:
                if minElements == 1:
                    ret = self + makeOptionalList(optElements)
                else:
                    ret = And([self] * minElements) + makeOptionalList(optElements)
            else:
                ret = makeOptionalList(optElements)
        else:
            if minElements == 1:
                ret = self
            else:
                ret = And([self] * minElements)
        return ret

    def __rmul__(self, other) -> \"ParserElement\":
        return self.__mul__(other)

    def __or__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``|`` operator - returns :class:`MatchFirst`
        \"\"\"
        if other is Ellipsis:
            return _PendingSkip(self, must_skip=True)

        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return MatchFirst([self, other])

    def __ror__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``|`` operator when left operand is not a :class:`ParserElement`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return other | self

    def __xor__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``^`` operator - returns :class:`Or`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return Or([self, other])

    def __rxor__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``^`` operator when left operand is not a :class:`ParserElement`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return other ^ self

    def __and__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``&`` operator - returns :class:`Each`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return Each([self, other])

    def __rand__(self, other) -> \"ParserElement\":
        \"\"\"
        Implementation of ``&`` operator when left operand is not a :class:`ParserElement`
        \"\"\"
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        if not isinstance(other, ParserElement):
            raise TypeError(
                \"Cannot combine element of type {} with ParserElement\".format(
                    type(other).__name__
                )
            )
        return other & self

    def __invert__(self) -> \"ParserElement\":
        \"\"\"
        Implementation of ``~`` operator - returns :class:`NotAny`
        \"\"\"
        return NotAny(self)

    # disable __iter__ to override legacy use of sequential access to __getitem__ to
    # iterate over a sequence
    __iter__ = None

    def __getitem__(self, key):
        \"\"\"
        use ``[]`` indexing notation as a short form for expression repetition:

        - ``expr[n]`` is equivalent to ``expr*n``
        - ``expr[m, n]`` is equivalent to ``expr*(m, n)``
        - ``expr[n, ...]`` or ``expr[n,]`` is equivalent
             to ``expr*n + ZeroOrMore(expr)``
             (read as \"at least n instances of ``expr``\")
        - ``expr[..., n]`` is equivalent to ``expr*(0, n)``
             (read as \"0 to n instances of ``expr``\")
        - ``expr[...]`` and ``expr[0, ...]`` are equivalent to ``ZeroOrMore(expr)``
        - ``expr[1, ...]`` is equivalent to ``OneOrMore(expr)``

        ``None`` may be used in place of ``...``.

        Note that ``expr[..., n]`` and ``expr[m, n]``do not raise an exception
        if more than ``n`` ``expr``s exist in the input stream.  If this behavior is
        desired, then write ``expr[..., n] + ~expr``.
        \"\"\"

        # convert single arg keys to tuples
        try:
            if isinstance(key, str_type):
                key = (key,)
            iter(key)
        except TypeError:
            key = (key, key)

        if len(key) > 2:
            raise TypeError(
                \"only 1 or 2 index arguments supported ({}{})\".format(
                    key[:5], \"... [{}]\".format(len(key)) if len(key) > 5 else \"\"
                )
            )

        # clip to 2 elements
        ret = self * tuple(key[:2])
        return ret

    def __call__(self, name: str = None) -> \"ParserElement\":
        \"\"\"
        Shortcut for :class:`set_results_name`, with ``list_all_matches=False``.

        If ``name`` is given with a trailing ``'*'`` character, then ``list_all_matches`` will be
        passed as ``True``.

        If ``name` is omitted, same as calling :class:`copy`.

        Example::

            # these are equivalent
            userdata = Word(alphas).set_results_name(\"name\") + Word(nums + \"-\").set_results_name(\"socsecno\")
            userdata = Word(alphas)(\"name\") + Word(nums + \"-\")(\"socsecno\")
        \"\"\"
        if name is not None:
            return self._setResultsName(name)
        else:
            return self.copy()

    def suppress(self) -> \"ParserElement\":
        \"\"\"
        Suppresses the output of this :class:`ParserElement`; useful to keep punctuation from
        cluttering up returned output.
        \"\"\"
        return Suppress(self)

    def ignore_whitespace(self, recursive: bool = True) -> \"ParserElement\":
        \"\"\"
        Enables the skipping of whitespace before matching the characters in the
        :class:`ParserElement`'s defined pattern.

        :param recursive: If ``True`` (the default), also enable whitespace skipping in child elements (if any)
        \"\"\"
        self.skipWhitespace = True
        return self

    def leave_whitespace(self, recursive: bool = True) -> \"ParserElement\":
        \"\"\"
        Disables the skipping of whitespace before matching the characters in the
        :class:`ParserElement`'s defined pattern.  This is normally only used internally by
        the pyparsing module, but may be needed in some whitespace-sensitive grammars.

        :param recursive: If true (the default), also disable whitespace skipping in child elements (if any)
        \"\"\"
        self.skipWhitespace = False
        return self

    def set_whitespace_chars(
        self, chars: Union[Set[str], str], copy_defaults: bool = False
    ) -> \"ParserElement\":
        \"\"\"
        Overrides the default whitespace chars
        \"\"\"
        self.skipWhitespace = True
        self.whiteChars = set(chars)
        self.copyDefaultWhiteChars = copy_defaults
        return self

    def parse_with_tabs(self) -> \"ParserElement\":
        \"\"\"
        Overrides default behavior to expand ``<TAB>`` s to spaces before parsing the input string.
        Must be called before ``parse_string`` when the input grammar contains elements that
        match ``<TAB>`` characters.
        \"\"\"
        self.keepTabs = True
        return self

    def ignore(self, other: \"ParserElement\") -> \"ParserElement\":
        \"\"\"
        Define expression to be ignored (e.g., comments) while doing pattern
        matching; may be called repeatedly, to define multiple comment or other
        ignorable patterns.

        Example::

            patt = Word(alphas)[1, ...]
            patt.parse_string('ablaj /* comment */ lskjd')
            # -> ['ablaj']

            patt.ignore(c_style_comment)
            patt.parse_string('ablaj /* comment */ lskjd')
            # -> ['ablaj', 'lskjd']
        \"\"\"
        import typing

        if isinstance(other, str_type):
            other = Suppress(other)

        if isinstance(other, Suppress):
            if other not in self.ignoreExprs:
                self.ignoreExprs.append(other)
        else:
            self.ignoreExprs.append(Suppress(other.copy()))
        return self

    def set_debug_actions(
        self,
        start_action: DebugStartAction,
        success_action: DebugSuccessAction,
        exception_action: DebugExceptionAction,
    ) -> \"ParserElement\":
        \"\"\"
        Customize display of debugging messages while doing pattern matching:

        - ``start_action`` - method to be called when an expression is about to be parsed;
          should have the signature ``fn(input_string: str, location: int, expression: ParserElement, cache_hit: bool)``

        - ``success_action`` - method to be called when an expression has successfully parsed;
          should have the signature ``fn(input_string: str, start_location: int, end_location: int, expression: ParserELement, parsed_tokens: ParseResults, cache_hit: bool)``

        - ``exception_action`` - method to be called when expression fails to parse;
          should have the signature ``fn(input_string: str, location: int, expression: ParserElement, exception: Exception, cache_hit: bool)``
        \"\"\"
        self.debugActions = self.DebugActions(
            start_action or _default_start_debug_action,
            success_action or _default_success_debug_action,
            exception_action or _default_exception_debug_action,
        )
        self.debug = True
        return self

    def set_debug(self, flag: bool = True) -> \"ParserElement\":
        \"\"\"
        Enable display of debugging messages while doing pattern matching.
        Set ``flag`` to ``True`` to enable, ``False`` to disable.

        Example::

            wd = Word(alphas).set_name(\"alphaword\")
            integer = Word(nums).set_name(\"numword\")
            term = wd | integer

            # turn on debugging for wd
            wd.set_debug()

            term[1, ...].parse_string(\"abc 123 xyz 890\")

        prints::

            Match alphaword at loc 0(1,1)
            Matched alphaword -> ['abc']
            Match alphaword at loc 3(1,4)
            Exception raised:Expected alphaword (at char 4), (line:1, col:5)
            Match alphaword at loc 7(1,8)
            Matched alphaword -> ['xyz']
            Match alphaword at loc 11(1,12)
            Exception raised:Expected alphaword (at char 12), (line:1, col:13)
            Match alphaword at loc 15(1,16)
            Exception raised:Expected alphaword (at char 15), (line:1, col:16)

        The output shown is that produced by the default debug actions - custom debug actions can be
        specified using :class:`set_debug_actions`. Prior to attempting
        to match the ``wd`` expression, the debugging message ``\"Match <exprname> at loc <n>(<line>,<col>)\"``
        is shown. Then if the parse succeeds, a ``\"Matched\"`` message is shown, or an ``\"Exception raised\"``
        message is shown. Also note the use of :class:`set_name` to assign a human-readable name to the expression,
        which makes debugging and exception messages easier to understand - for instance, the default
        name created for the :class:`Word` expression without calling ``set_name`` is ``\"W:(A-Za-z)\"``.
        \"\"\"
        if flag:
            self.set_debug_actions(
                _default_start_debug_action,
                _default_success_debug_action,
                _default_exception_debug_action,
            )
        else:
            self.debug = False
        return self

    @property
    def default_name(self) -> str:
        if self._defaultName is None:
            self._defaultName = self._generateDefaultName()
        return self._defaultName

    @abstractmethod
    def _generateDefaultName(self):
        \"\"\"
        Child classes must define this method, which defines how the ``default_name`` is set.
        \"\"\"

    def set_name(self, name: str) -> \"ParserElement\":
        \"\"\"
        Define name for this expression, makes debugging and exception messages clearer.
        Example::
            Word(nums).parse_string(\"ABC\")  # -> Exception: Expected W:(0-9) (at char 0), (line:1, col:1)
            Word(nums).set_name(\"integer\").parse_string(\"ABC\")  # -> Exception: Expected integer (at char 0), (line:1, col:1)
        \"\"\"
        self.customName = name
        self.errmsg = \"Expected \" + self.name
        if __diag__.enable_debug_on_named_expressions:
            self.set_debug()
        return self

    @property
    def name(self) -> str:
        # This will use a user-defined name if available, but otherwise defaults back to the auto-generated name
        return self.customName if self.customName is not None else self.default_name

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return str(self)

    def streamline(self) -> \"ParserElement\":
        self.streamlined = True
        self._defaultName = None
        return self

    def recurse(self) -> Sequence[\"ParserElement\"]:
        return []

    def _checkRecursion(self, parseElementList):
        subRecCheckList = parseElementList[:] + [self]
        for e in self.recurse():
            e._checkRecursion(subRecCheckList)

    def validate(self, validateTrace=None) -> None:
        \"\"\"
        Check defined expressions for valid structure, check for infinite recursive definitions.
        \"\"\"
        self._checkRecursion([])

    def parse_file(
        self,
        file_or_filename: Union[str, Path, TextIO],
        encoding: str = \"utf-8\",
        parse_all: bool = False,
        *,
        parseAll: bool = False,
    ) -> ParseResults:
        \"\"\"
        Execute the parse expression on the given file or filename.
        If a filename is specified (instead of a file object),
        the entire file is opened, read, and closed before parsing.
        \"\"\"
        parseAll = parseAll or parse_all
        try:
            file_contents = file_or_filename.read()
        except AttributeError:
            with open(file_or_filename, \"r\", encoding=encoding) as f:
                file_contents = f.read()
        try:
            return self.parse_string(file_contents, parseAll)
        except ParseBaseException as exc:
            if ParserElement.verbose_stacktrace:
                raise
            else:
                # catch and re-raise exception from here, clears out pyparsing internal stack trace
                raise exc.with_traceback(None)

    def __eq__(self, other):
        if self is other:
            return True
        elif isinstance(other, str_type):
            return self.matches(other, parse_all=True)
        elif isinstance(other, ParserElement):
            return vars(self) == vars(other)
        return False

    def __hash__(self):
        return id(self)

    def matches(
        self, test_string: str, parse_all: bool = True, *, parseAll: bool = True
    ) -> bool:
        \"\"\"
        Method for quick testing of a parser against a test string. Good for simple
        inline microtests of sub expressions while building up larger parser.

        Parameters:
        - ``test_string`` - to test against this expression for a match
        - ``parse_all`` - (default= ``True``) - flag to pass to :class:`parse_string` when running tests

        Example::

            expr = Word(nums)
            assert expr.matches(\"100\")
        \"\"\"
        parseAll = parseAll and parse_all
        try:
            self.parse_string(str(test_string), parse_all=parseAll)
            return True
        except ParseBaseException:
            return False

    def run_tests(
        self,
        tests: Union[str, List[str]],
        parse_all: bool = True,
        comment: typing.Optional[Union[\"ParserElement\", str]] = \"#\",
        full_dump: bool = True,
        print_results: bool = True,
        failure_tests: bool = False,
        post_parse: Callable[[str, ParseResults], str] = None,
        file: typing.Optional[TextIO] = None,
        with_line_numbers: bool = False,
        *,
        parseAll: bool = True,
        fullDump: bool = True,
        printResults: bool = True,
        failureTests: bool = False,
        postParse: Callable[[str, ParseResults], str] = None,
    ) -> Tuple[bool, List[Tuple[str, Union[ParseResults, Exception]]]]:
        \"\"\"
        Execute the parse expression on a series of test strings, showing each
        test, the parsed results or where the parse failed. Quick and easy way to
        run a parse expression against a list of sample strings.

        Parameters:
        - ``tests`` - a list of separate test strings, or a multiline string of test strings
        - ``parse_all`` - (default= ``True``) - flag to pass to :class:`parse_string` when running tests
        - ``comment`` - (default= ``'#'``) - expression for indicating embedded comments in the test
          string; pass None to disable comment filtering
        - ``full_dump`` - (default= ``True``) - dump results as list followed by results names in nested outline;
          if False, only dump nested list
        - ``print_results`` - (default= ``True``) prints test output to stdout
        - ``failure_tests`` - (default= ``False``) indicates if these tests are expected to fail parsing
        - ``post_parse`` - (default= ``None``) optional callback for successful parse results; called as
          `fn(test_string, parse_results)` and returns a string to be added to the test output
        - ``file`` - (default= ``None``) optional file-like object to which test output will be written;
          if None, will default to ``sys.stdout``
        - ``with_line_numbers`` - default= ``False``) show test strings with line and column numbers

        Returns: a (success, results) tuple, where success indicates that all tests succeeded
        (or failed if ``failure_tests`` is True), and the results contain a list of lines of each
        test's output

        Example::

            number_expr = pyparsing_common.number.copy()

            result = number_expr.run_tests('''
                # unsigned integer
                100
                # negative integer
                -100
                # float with scientific notation
                6.02e23
                # integer with scientific notation
                1e-12
                ''')
            print(\"Success\" if result[0] else \"Failed!\")

            result = number_expr.run_tests('''
                # stray character
                100Z
                # missing leading digit before '.'
                -.100
                # too many '.'
                3.14.159
                ''', failure_tests=True)
            print(\"Success\" if result[0] else \"Failed!\")

        prints::

            # unsigned integer
            100
            [100]

            # negative integer
            -100
            [-100]

            # float with scientific notation
            6.02e23
            [6.02e+23]

            # integer with scientific notation
            1e-12
            [1e-12]

            Success

            # stray character
            100Z
               ^
            FAIL: Expected end of text (at char 3), (line:1, col:4)

            # missing leading digit before '.'
            -.100
            ^
            FAIL: Expected {real number with scientific notation | real number | signed integer} (at char 0), (line:1, col:1)

            # too many '.'
            3.14.159
                ^
            FAIL: Expected end of text (at char 4), (line:1, col:5)

            Success

        Each test string must be on a single line. If you want to test a string that spans multiple
        lines, create a test like this::

            expr.run_tests(r\"this is a test\\\\n of strings that spans \\\\n 3 lines\")

        (Note that this is a raw string literal, you must include the leading ``'r'``.)
        \"\"\"
        from .testing import pyparsing_test

        parseAll = parseAll and parse_all
        fullDump = fullDump and full_dump
        printResults = printResults and print_results
        failureTests = failureTests or failure_tests
        postParse = postParse or post_parse
        if isinstance(tests, str_type):
            line_strip = type(tests).strip
            tests = [line_strip(test_line) for test_line in tests.rstrip().splitlines()]
        if isinstance(comment, str_type):
            comment = Literal(comment)
        if file is None:
            file = sys.stdout
        print_ = file.write

        result: Union[ParseResults, Exception]
        allResults = []
        comments = []
        success = True
        NL = Literal(r\"\\n\").add_parse_action(replace_with(\"\\n\")).ignore(quoted_string)
        BOM = \"\\ufeff\"
        for t in tests:
            if comment is not None and comment.matches(t, False) or comments and not t:
                comments.append(
                    pyparsing_test.with_line_numbers(t) if with_line_numbers else t
                )
                continue
            if not t:
                continue
            out = [
                \"\\n\" + \"\\n\".join(comments) if comments else \"\",
                pyparsing_test.with_line_numbers(t) if with_line_numbers else t,
            ]
            comments = []
            try:
                # convert newline marks to actual newlines, and strip leading BOM if present
                t = NL.transform_string(t.lstrip(BOM))
                result = self.parse_string(t, parse_all=parseAll)
            except ParseBaseException as pe:
                fatal = \"(FATAL)\" if isinstance(pe, ParseFatalException) else \"\"
                out.append(pe.explain())
                out.append(\"FAIL: \" + str(pe))
                if ParserElement.verbose_stacktrace:
                    out.extend(traceback.format_tb(pe.__traceback__))
                success = success and failureTests
                result = pe
            except Exception as exc:
                out.append(\"FAIL-EXCEPTION: {}: {}\".format(type(exc).__name__, exc))
                if ParserElement.verbose_stacktrace:
                    out.extend(traceback.format_tb(exc.__traceback__))
                success = success and failureTests
                result = exc
            else:
                success = success and not failureTests
                if postParse is not None:
                    try:
                        pp_value = postParse(t, result)
                        if pp_value is not None:
                            if isinstance(pp_value, ParseResults):
                                out.append(pp_value.dump())
                            else:
                                out.append(str(pp_value))
                        else:
                            out.append(result.dump())
                    except Exception as e:
                        out.append(result.dump(full=fullDump))
                        out.append(
                            \"{} failed: {}: {}\".format(
                                postParse.__name__, type(e).__name__, e
                            )
                        )
                else:
                    out.append(result.dump(full=fullDump))
            out.append(\"\")

            if printResults:
                print_(\"\\n\".join(out))

            allResults.append((t, result))

        return success, allResults

    def create_diagram(
        self,
        output_html: Union[TextIO, Path, str],
        vertical: int = 3,
        show_results_names: bool = False,
        show_groups: bool = False,
        **kwargs,
    ) -> None:
        \"\"\"
        Create a railroad diagram for the parser.

        Parameters:
        - output_html (str or file-like object) - output target for generated
          diagram HTML
        - vertical (int) - threshold for formatting multiple alternatives vertically
          instead of horizontally (default=3)
        - show_results_names - bool flag whether diagram should show annotations for
          defined results names
        - show_groups - bool flag whether groups should be highlighted with an unlabeled surrounding box
        Additional diagram-formatting keyword arguments can also be included;
        see railroad.Diagram class.
        \"\"\"

        try:
            from .diagram import to_railroad, railroad_to_html
        except ImportError as ie:
            raise Exception(
                \"must ``pip install pyparsing[diagrams]`` to generate parser railroad diagrams\"
            ) from ie

        self.streamline()

        railroad = to_railroad(
            self,
            vertical=vertical,
            show_results_names=show_results_names,
            show_groups=show_groups,
            diagram_kwargs=kwargs,
        )
        if isinstance(output_html, (str, Path)):
            with open(output_html, \"w\", encoding=\"utf-8\") as diag_file:
                diag_file.write(railroad_to_html(railroad))
        else:
            # we were passed a file-like object, just write to it
            output_html.write(railroad_to_html(railroad))

    setDefaultWhitespaceChars = set_default_whitespace_chars
    inlineLiteralsUsing = inline_literals_using
    setResultsName = set_results_name
    setBreak = set_break
    setParseAction = set_parse_action
    addParseAction = add_parse_action
    addCondition = add_condition
    setFailAction = set_fail_action
    tryParse = try_parse
    canParseNext = can_parse_next
    resetCache = reset_cache
    enableLeftRecursion = enable_left_recursion
    enablePackrat = enable_packrat
    parseString = parse_string
    scanString = scan_string
    searchString = search_string
    transformString = transform_string
    setWhitespaceChars = set_whitespace_chars
    parseWithTabs = parse_with_tabs
    setDebugActions = set_debug_actions
    setDebug = set_debug
    defaultName = default_name
    setName = set_name
    parseFile = parse_file
    runTests = run_tests
    ignoreWhitespace = ignore_whitespace
    leaveWhitespace = leave_whitespace


class _PendingSkip(ParserElement):
    # internal placeholder class to hold a place were '...' is added to a parser element,
    # once another ParserElement is added, this placeholder will be replaced with a SkipTo
    def __init__(self, expr: ParserElement, must_skip: bool = False):
        super().__init__()
        self.anchor = expr
        self.must_skip = must_skip

    def _generateDefaultName(self):
        return str(self.anchor + Empty()).replace(\"Empty\", \"...\")

    def __add__(self, other) -> \"ParserElement\":
        skipper = SkipTo(other).set_name(\"...\")(\"_skipped*\")
        if self.must_skip:

            def must_skip(t):
                if not t._skipped or t._skipped.as_list() == [\"\"]:
                    del t[0]
                    t.pop(\"_skipped\", None)

            def show_skip(t):
                if t._skipped.as_list()[-1:] == [\"\"]:
                    t.pop(\"_skipped\")
                    t[\"_skipped\"] = \"missing <\" + repr(self.anchor) + \">\"

            return (
                self.anchor + skipper().add_parse_action(must_skip)
                | skipper().add_parse_action(show_skip)
            ) + other

        return self.anchor + skipper + other

    def __repr__(self):
        return self.defaultName

    def parseImpl(self, *args):
        raise Exception(
            \"use of `...` expression without following SkipTo target expression\"
        )


class Token(ParserElement):
    \"\"\"Abstract :class:`ParserElement` subclass, for defining atomic
    matching patterns.
    \"\"\"

    def __init__(self):
        super().__init__(savelist=False)

    def _generateDefaultName(self):
        return type(self).__name__


class Empty(Token):
    \"\"\"
    An empty token, will always match.
    \"\"\"

    def __init__(self):
        super().__init__()
        self.mayReturnEmpty = True
        self.mayIndexError = False


class NoMatch(Token):
    \"\"\"
    A token that will never match.
    \"\"\"

    def __init__(self):
        super().__init__()
        self.mayReturnEmpty = True
        self.mayIndexError = False
        self.errmsg = \"Unmatchable token\"

    def parseImpl(self, instring, loc, doActions=True):
        raise ParseException(instring, loc, self.errmsg, self)


class Literal(Token):
    \"\"\"
    Token to exactly match a specified string.

    Example::

        Literal('blah').parse_string('blah')  # -> ['blah']
        Literal('blah').parse_string('blahfooblah')  # -> ['blah']
        Literal('blah').parse_string('bla')  # -> Exception: Expected \"blah\"

    For case-insensitive matching, use :class:`CaselessLiteral`.

    For keyword matching (force word break before and after the matched string),
    use :class:`Keyword` or :class:`CaselessKeyword`.
    \"\"\"

    def __init__(self, match_string: str = \"\", *, matchString: str = \"\"):
        super().__init__()
        match_string = matchString or match_string
        self.match = match_string
        self.matchLen = len(match_string)
        try:
            self.firstMatchChar = match_string[0]
        except IndexError:
            raise ValueError(\"null string passed to Literal; use Empty() instead\")
        self.errmsg = \"Expected \" + self.name
        self.mayReturnEmpty = False
        self.mayIndexError = False

        # Performance tuning: modify __class__ to select
        # a parseImpl optimized for single-character check
        if self.matchLen == 1 and type(self) is Literal:
            self.__class__ = _SingleCharLiteral

    def _generateDefaultName(self):
        return repr(self.match)

    def parseImpl(self, instring, loc, doActions=True):
        if instring[loc] == self.firstMatchChar and instring.startswith(
            self.match, loc
        ):
            return loc + self.matchLen, self.match
        raise ParseException(instring, loc, self.errmsg, self)


class _SingleCharLiteral(Literal):
    def parseImpl(self, instring, loc, doActions=True):
        if instring[loc] == self.firstMatchChar:
            return loc + 1, self.match
        raise ParseException(instring, loc, self.errmsg, self)


ParserElement._literalStringClass = Literal


class Keyword(Token):
    \"\"\"
    Token to exactly match a specified string as a keyword, that is,
    it must be immediately followed by a non-keyword character.  Compare
    with :class:`Literal`:

    - ``Literal(\"if\")`` will match the leading ``'if'`` in
      ``'ifAndOnlyIf'``.
    - ``Keyword(\"if\")`` will not; it will only match the leading
      ``'if'`` in ``'if x=1'``, or ``'if(y==2)'``

    Accepts two optional constructor arguments in addition to the
    keyword string:

    - ``identChars`` is a string of characters that would be valid
      identifier characters, defaulting to all alphanumerics + \"_\" and
      \"$\"
    - ``caseless`` allows case-insensitive matching, default is ``False``.

    Example::

        Keyword(\"start\").parse_string(\"start\")  # -> ['start']
        Keyword(\"start\").parse_string(\"starting\")  # -> Exception

    For case-insensitive matching, use :class:`CaselessKeyword`.
    \"\"\"

    DEFAULT_KEYWORD_CHARS = alphanums + \"_$\"

    def __init__(
        self,
        match_string: str = \"\",
        ident_chars: typing.Optional[str] = None,
        caseless: bool = False,
        *,
        matchString: str = \"\",
        identChars: typing.Optional[str] = None,
    ):
        super().__init__()
        identChars = identChars or ident_chars
        if identChars is None:
            identChars = Keyword.DEFAULT_KEYWORD_CHARS
        match_string = matchString or match_string
        self.match = match_string
        self.matchLen = len(match_string)
        try:
            self.firstMatchChar = match_string[0]
        except IndexError:
            raise ValueError(\"null string passed to Keyword; use Empty() instead\")
        self.errmsg = \"Expected {} {}\".format(type(self).__name__, self.name)
        self.mayReturnEmpty = False
        self.mayIndexError = False
        self.caseless = caseless
        if caseless:
            self.caselessmatch = match_string.upper()
            identChars = identChars.upper()
        self.identChars = set(identChars)

    def _generateDefaultName(self):
        return repr(self.match)

    def parseImpl(self, instring, loc, doActions=True):
        errmsg = self.errmsg
        errloc = loc
        if self.caseless:
            if instring[loc : loc + self.matchLen].upper() == self.caselessmatch:
                if loc == 0 or instring[loc - 1].upper() not in self.identChars:
                    if (
                        loc >= len(instring) - self.matchLen
                        or instring[loc + self.matchLen].upper() not in self.identChars
                    ):
                        return loc + self.matchLen, self.match
                    else:
                        # followed by keyword char
                        errmsg += \", was immediately followed by keyword character\"
                        errloc = loc + self.matchLen
                else:
                    # preceded by keyword char
                    errmsg += \", keyword was immediately preceded by keyword character\"
                    errloc = loc - 1
            # else no match just raise plain exception

        else:
            if (
                instring[loc] == self.firstMatchChar
                and self.matchLen == 1
                or instring.startswith(self.match, loc)
            ):
                if loc == 0 or instring[loc - 1] not in self.identChars:
                    if (
                        loc >= len(instring) - self.matchLen
                        or instring[loc + self.matchLen] not in self.identChars
                    ):
                        return loc + self.matchLen, self.match
                    else:
                        # followed by keyword char
                        errmsg += (
                            \", keyword was immediately followed by keyword character\"
                        )
                        errloc = loc + self.matchLen
                else:
                    # preceded by keyword char
                    errmsg += \", keyword was immediately preceded by keyword character\"
                    errloc = loc - 1
            # else no match just raise plain exception

        raise ParseException(instring, errloc, errmsg, self)

    @staticmethod
    def set_default_keyword_chars(chars) -> None:
        \"\"\"
        Overrides the default characters used by :class:`Keyword` expressions.
        \"\"\"
        Keyword.DEFAULT_KEYWORD_CHARS = chars

    setDefaultKeywordChars = set_default_keyword_chars


class CaselessLiteral(Literal):
    \"\"\"
    Token to match a specified string, ignoring case of letters.
    Note: the matched results will always be in the case of the given
    match string, NOT the case of the input text.

    Example::

        CaselessLiteral(\"CMD\")[1, ...].parse_string(\"cmd CMD Cmd10\")
        # -> ['CMD', 'CMD', 'CMD']

    (Contrast with example for :class:`CaselessKeyword`.)
    \"\"\"

    def __init__(self, match_string: str = \"\", *, matchString: str = \"\"):
        match_string = matchString or match_string
        super().__init__(match_string.upper())
        # Preserve the defining literal.
        self.returnString = match_string
        self.errmsg = \"Expected \" + self.name

    def parseImpl(self, instring, loc, doActions=True):
        if instring[loc : loc + self.matchLen].upper() == self.match:
            return loc + self.matchLen, self.returnString
        raise ParseException(instring, loc, self.errmsg, self)


class CaselessKeyword(Keyword):
    \"\"\"
    Caseless version of :class:`Keyword`.

    Example::

        CaselessKeyword(\"CMD\")[1, ...].parse_string(\"cmd CMD Cmd10\")
        # -> ['CMD', 'CMD']

    (Contrast with example for :class:`CaselessLiteral`.)
    \"\"\"

    def __init__(
        self,
        match_string: str = \"\",
        ident_chars: typing.Optional[str] = None,
        *,
        matchString: str = \"\",
        identChars: typing.Optional[str] = None,
    ):
        identChars = identChars or ident_chars
        match_string = matchString or match_string
        super().__init__(match_string, identChars, caseless=True)


class CloseMatch(Token):
    \"\"\"A variation on :class:`Literal` which matches \"close\" matches,
    that is, strings with at most 'n' mismatching characters.
    :class:`CloseMatch` takes parameters:

    - ``match_string`` - string to be matched
    - ``caseless`` - a boolean indicating whether to ignore casing when comparing characters
    - ``max_mismatches`` - (``default=1``) maximum number of
      mismatches allowed to count as a match

    The results from a successful parse will contain the matched text
    from the input string and the following named results:

    - ``mismatches`` - a list of the positions within the
      match_string where mismatches were found
    - ``original`` - the original match_string used to compare
      against the input string

    If ``mismatches`` is an empty list, then the match was an exact
    match.

    Example::

        patt = CloseMatch(\"ATCATCGAATGGA\")
        patt.parse_string(\"ATCATCGAAXGGA\") # -> (['ATCATCGAAXGGA'], {'mismatches': [[9]], 'original': ['ATCATCGAATGGA']})
        patt.parse_string(\"ATCAXCGAAXGGA\") # -> Exception: Expected 'ATCATCGAATGGA' (with up to 1 mismatches) (at char 0), (line:1, col:1)

        # exact match
        patt.parse_string(\"ATCATCGAATGGA\") # -> (['ATCATCGAATGGA'], {'mismatches': [[]], 'original': ['ATCATCGAATGGA']})

        # close match allowing up to 2 mismatches
        patt = CloseMatch(\"ATCATCGAATGGA\", max_mismatches=2)
        patt.parse_string(\"ATCAXCGAAXGGA\") # -> (['ATCAXCGAAXGGA'], {'mismatches': [[4, 9]], 'original': ['ATCATCGAATGGA']})
    \"\"\"

    def __init__(
        self,
        match_string: str,
        max_mismatches: int = None,
        *,
        maxMismatches: int = 1,
        caseless=False,
    ):
        maxMismatches = max_mismatches if max_mismatches is not None else maxMismatches
        super().__init__()
        self.match_string = match_string
        self.maxMismatches = maxMismatches
        self.errmsg = \"Expected {!r} (with up to {} mismatches)\".format(
            self.match_string, self.maxMismatches
        )
        self.caseless = caseless
        self.mayIndexError = False
        self.mayReturnEmpty = False

    def _generateDefaultName(self):
        return \"{}:{!r}\".format(type(self).__name__, self.match_string)

    def parseImpl(self, instring, loc, doActions=True):
        start = loc
        instrlen = len(instring)
        maxloc = start + len(self.match_string)

        if maxloc <= instrlen:
            match_string = self.match_string
            match_stringloc = 0
            mismatches = []
            maxMismatches = self.maxMismatches

            for match_stringloc, s_m in enumerate(
                zip(instring[loc:maxloc], match_string)
            ):
                src, mat = s_m
                if self.caseless:
                    src, mat = src.lower(), mat.lower()

                if src != mat:
                    mismatches.append(match_stringloc)
                    if len(mismatches) > maxMismatches:
                        break
            else:
                loc = start + match_stringloc + 1
                results = ParseResults([instring[start:loc]])
                results[\"original\"] = match_string
                results[\"mismatches\"] = mismatches
                return loc, results

        raise ParseException(instring, loc, self.errmsg, self)


class Word(Token):
    \"\"\"Token for matching words composed of allowed character sets.
    Parameters:
    - ``init_chars`` - string of all characters that should be used to
      match as a word; \"ABC\" will match \"AAA\", \"ABAB\", \"CBAC\", etc.;
      if ``body_chars`` is also specified, then this is the string of
      initial characters
    - ``body_chars`` - string of characters that
      can be used for matching after a matched initial character as
      given in ``init_chars``; if omitted, same as the initial characters
      (default=``None``)
    - ``min`` - minimum number of characters to match (default=1)
    - ``max`` - maximum number of characters to match (default=0)
    - ``exact`` - exact number of characters to match (default=0)
    - ``as_keyword`` - match as a keyword (default=``False``)
    - ``exclude_chars`` - characters that might be
      found in the input ``body_chars`` string but which should not be
      accepted for matching ;useful to define a word of all
      printables except for one or two characters, for instance
      (default=``None``)

    :class:`srange` is useful for defining custom character set strings
    for defining :class:`Word` expressions, using range notation from
    regular expression character sets.

    A common mistake is to use :class:`Word` to match a specific literal
    string, as in ``Word(\"Address\")``. Remember that :class:`Word`
    uses the string argument to define *sets* of matchable characters.
    This expression would match \"Add\", \"AAA\", \"dAred\", or any other word
    made up of the characters 'A', 'd', 'r', 'e', and 's'. To match an
    exact literal string, use :class:`Literal` or :class:`Keyword`.

    pyparsing includes helper strings for building Words:

    - :class:`alphas`
    - :class:`nums`
    - :class:`alphanums`
    - :class:`hexnums`
    - :class:`alphas8bit` (alphabetic characters in ASCII range 128-255
      - accented, tilded, umlauted, etc.)
    - :class:`punc8bit` (non-alphabetic characters in ASCII range
      128-255 - currency, symbols, superscripts, diacriticals, etc.)
    - :class:`printables` (any non-whitespace character)

    ``alphas``, ``nums``, and ``printables`` are also defined in several
    Unicode sets - see :class:`pyparsing_unicode``.

    Example::

        # a word composed of digits
        integer = Word(nums) # equivalent to Word(\"0123456789\") or Word(srange(\"0-9\"))

        # a word with a leading capital, and zero or more lowercase
        capital_word = Word(alphas.upper(), alphas.lower())

        # hostnames are alphanumeric, with leading alpha, and '-'
        hostname = Word(alphas, alphanums + '-')

        # roman numeral (not a strict parser, accepts invalid mix of characters)
        roman = Word(\"IVXLCDM\")

        # any string of non-whitespace characters, except for ','
        csv_value = Word(printables, exclude_chars=\",\")
    \"\"\"

    def __init__(
        self,
        init_chars: str = \"\",
        body_chars: typing.Optional[str] = None,
        min: int = 1,
        max: int = 0,
        exact: int = 0,
        as_keyword: bool = False,
        exclude_chars: typing.Optional[str] = None,
        *,
        initChars: typing.Optional[str] = None,
        bodyChars: typing.Optional[str] = None,
        asKeyword: bool = False,
        excludeChars: typing.Optional[str] = None,
    ):
        initChars = initChars or init_chars
        bodyChars = bodyChars or body_chars
        asKeyword = asKeyword or as_keyword
        excludeChars = excludeChars or exclude_chars
        super().__init__()
        if not initChars:
            raise ValueError(
                \"invalid {}, initChars cannot be empty string\".format(
                    type(self).__name__
                )
            )

        initChars = set(initChars)
        self.initChars = initChars
        if excludeChars:
            excludeChars = set(excludeChars)
            initChars -= excludeChars
            if bodyChars:
                bodyChars = set(bodyChars) - excludeChars
        self.initCharsOrig = \"\".join(sorted(initChars))

        if bodyChars:
            self.bodyCharsOrig = \"\".join(sorted(bodyChars))
            self.bodyChars = set(bodyChars)
        else:
            self.bodyCharsOrig = \"\".join(sorted(initChars))
            self.bodyChars = set(initChars)

        self.maxSpecified = max > 0

        if min < 1:
            raise ValueError(
                \"cannot specify a minimum length < 1; use Opt(Word()) if zero-length word is permitted\"
            )

        self.minLen = min

        if max > 0:
            self.maxLen = max
        else:
            self.maxLen = _MAX_INT

        if exact > 0:
            self.maxLen = exact
            self.minLen = exact

        self.errmsg = \"Expected \" + self.name
        self.mayIndexError = False
        self.asKeyword = asKeyword

        # see if we can make a regex for this Word
        if \" \" not in self.initChars | self.bodyChars and (min == 1 and exact == 0):
            if self.bodyChars == self.initChars:
                if max == 0:
                    repeat = \"+\"
                elif max == 1:
                    repeat = \"\"
                else:
                    repeat = \"{{{},{}}}\".format(
                        self.minLen, \"\" if self.maxLen == _MAX_INT else self.maxLen
                    )
                self.reString = \"[{}]{}\".format(
                    _collapse_string_to_ranges(self.initChars),
                    repeat,
                )
            elif len(self.initChars) == 1:
                if max == 0:
                    repeat = \"*\"
                else:
                    repeat = \"{{0,{}}}\".format(max - 1)
                self.reString = \"{}[{}]{}\".format(
                    re.escape(self.initCharsOrig),
                    _collapse_string_to_ranges(self.bodyChars),
                    repeat,
                )
            else:
                if max == 0:
                    repeat = \"*\"
                elif max == 2:
                    repeat = \"\"
                else:
                    repeat = \"{{0,{}}}\".format(max - 1)
                self.reString = \"[{}][{}]{}\".format(
                    _collapse_string_to_ranges(self.initChars),
                    _collapse_string_to_ranges(self.bodyChars),
                    repeat,
                )
            if self.asKeyword:
                self.reString = r\"\\b\" + self.reString + r\"\\b\"

            try:
                self.re = re.compile(self.reString)
            except re.error:
                self.re = None
            else:
                self.re_match = self.re.match
                self.__class__ = _WordRegex

    def _generateDefaultName(self):
        def charsAsStr(s):
            max_repr_len = 16
            s = _collapse_string_to_ranges(s, re_escape=False)
            if len(s) > max_repr_len:
                return s[: max_repr_len - 3] + \"...\"
            else:
                return s

        if self.initChars != self.bodyChars:
            base = \"W:({}, {})\".format(
                charsAsStr(self.initChars), charsAsStr(self.bodyChars)
            )
        else:
            base = \"W:({})\".format(charsAsStr(self.initChars))

        # add length specification
        if self.minLen > 1 or self.maxLen != _MAX_INT:
            if self.minLen == self.maxLen:
                if self.minLen == 1:
                    return base[2:]
                else:
                    return base + \"{{{}}}\".format(self.minLen)
            elif self.maxLen == _MAX_INT:
                return base + \"{{{},...}}\".format(self.minLen)
            else:
                return base + \"{{{},{}}}\".format(self.minLen, self.maxLen)
        return base

    def parseImpl(self, instring, loc, doActions=True):
        if instring[loc] not in self.initChars:
            raise ParseException(instring, loc, self.errmsg, self)

        start = loc
        loc += 1
        instrlen = len(instring)
        bodychars = self.bodyChars
        maxloc = start + self.maxLen
        maxloc = min(maxloc, instrlen)
        while loc < maxloc and instring[loc] in bodychars:
            loc += 1

        throwException = False
        if loc - start < self.minLen:
            throwException = True
        elif self.maxSpecified and loc < instrlen and instring[loc] in bodychars:
            throwException = True
        elif self.asKeyword:
            if (
                start > 0
                and instring[start - 1] in bodychars
                or loc < instrlen
                and instring[loc] in bodychars
            ):
                throwException = True

        if throwException:
            raise ParseException(instring, loc, self.errmsg, self)

        return loc, instring[start:loc]


class _WordRegex(Word):
    def parseImpl(self, instring, loc, doActions=True):
        result = self.re_match(instring, loc)
        if not result:
            raise ParseException(instring, loc, self.errmsg, self)

        loc = result.end()
        return loc, result.group()


class Char(_WordRegex):
    \"\"\"A short-cut class for defining :class:`Word` ``(characters, exact=1)``,
    when defining a match of any single character in a string of
    characters.
    \"\"\"

    def __init__(
        self,
        charset: str,
        as_keyword: bool = False,
        exclude_chars: typing.Optional[str] = None,
        *,
        asKeyword: bool = False,
        excludeChars: typing.Optional[str] = None,
    ):
        asKeyword = asKeyword or as_keyword
        excludeChars = excludeChars or exclude_chars
        super().__init__(
            charset, exact=1, asKeyword=asKeyword, excludeChars=excludeChars
        )
        self.reString = \"[{}]\".format(_collapse_string_to_ranges(self.initChars))
        if asKeyword:
            self.reString = r\"\\b{}\\b\".format(self.reString)
        self.re = re.compile(self.reString)
        self.re_match = self.re.match


class Regex(Token):
    r\"\"\"Token for matching strings that match a given regular
    expression. Defined with string specifying the regular expression in
    a form recognized by the stdlib Python  `re module <https://docs.python.org/3/library/re.html>`_.
    If the given regex contains named groups (defined using ``(?P<name>...)``),
    these will be preserved as named :class:`ParseResults`.

    If instead of the Python stdlib ``re`` module you wish to use a different RE module
    (such as the ``regex`` module), you can do so by building your ``Regex`` object with
    a compiled RE that was compiled using ``regex``.

    Example::

        realnum = Regex(r\"[+-]?\\d+\\.\\d*\")
        # ref: https://stackoverflow.com/questions/267399/how-do-you-match-only-valid-roman-numerals-with-a-regular-expression
        roman = Regex(r\"M{0,4}(CM|CD|D?{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})\")

        # named fields in a regex will be returned as named results
        date = Regex(r'(?P<year>\\d{4})-(?P<month>\\d\\d?)-(?P<day>\\d\\d?)')

        # the Regex class will accept re's compiled using the regex module
        import regex
        parser = pp.Regex(regex.compile(r'[0-9]'))
    \"\"\"

    def __init__(
        self,
        pattern: Any,
        flags: Union[re.RegexFlag, int] = 0,
        as_group_list: bool = False,
        as_match: bool = False,
        *,
        asGroupList: bool = False,
        asMatch: bool = False,
    ):
        \"\"\"The parameters ``pattern`` and ``flags`` are passed
        to the ``re.compile()`` function as-is. See the Python
        `re module <https://docs.python.org/3/library/re.html>`_ module for an
        explanation of the acceptable patterns and flags.
        \"\"\"
        super().__init__()
        asGroupList = asGroupList or as_group_list
        asMatch = asMatch or as_match

        if isinstance(pattern, str_type):
            if not pattern:
                raise ValueError(\"null string passed to Regex; use Empty() instead\")

            self._re = None
            self.reString = self.pattern = pattern
            self.flags = flags

        elif hasattr(pattern, \"pattern\") and hasattr(pattern, \"match\"):
            self._re = pattern
            self.pattern = self.reString = pattern.pattern
            self.flags = flags

        else:
            raise TypeError(
                \"Regex may only be constructed with a string or a compiled RE object\"
            )

        self.errmsg = \"Expected \" + self.name
        self.mayIndexError = False
        self.asGroupList = asGroupList
        self.asMatch = asMatch
        if self.asGroupList:
            self.parseImpl = self.parseImplAsGroupList
        if self.asMatch:
            self.parseImpl = self.parseImplAsMatch

    @cached_property
    def re(self):
        if self._re:
            return self._re
        else:
            try:
                return re.compile(self.pattern, self.flags)
            except re.error:
                raise ValueError(
                    \"invalid pattern ({!r}) passed to Regex\".format(self.pattern)
                )

    @cached_property
    def re_match(self):
        return self.re.match

    @cached_property
    def mayReturnEmpty(self):
        return self.re_match(\"\") is not None

    def _generateDefaultName(self):
        return \"Re:({})\".format(repr(self.pattern).replace(\"\\\\\\\\\", \"\\\\\"))

    def parseImpl(self, instring, loc, doActions=True):
        result = self.re_match(instring, loc)
        if not result:
            raise ParseException(instring, loc, self.errmsg, self)

        loc = result.end()
        ret = ParseResults(result.group())
        d = result.groupdict()
        if d:
            for k, v in d.items():
                ret[k] = v
        return loc, ret

    def parseImplAsGroupList(self, instring, loc, doActions=True):
        result = self.re_match(instring, loc)
        if not result:
            raise ParseException(instring, loc, self.errmsg, self)

        loc = result.end()
        ret = result.groups()
        return loc, ret

    def parseImplAsMatch(self, instring, loc, doActions=True):
        result = self.re_match(instring, loc)
        if not result:
            raise ParseException(instring, loc, self.errmsg, self)

        loc = result.end()
        ret = result
        return loc, ret

    def sub(self, repl: str) -> ParserElement:
        r\"\"\"
        Return :class:`Regex` with an attached parse action to transform the parsed
        result as if called using `re.sub(expr, repl, string) <https://docs.python.org/3/library/re.html#re.sub>`_.

        Example::

            make_html = Regex(r\"(\\w+):(.*?):\").sub(r\"<\\1>\\2</\\1>\")
            print(make_html.transform_string(\"h1:main title:\"))
            # prints \"<h1>main title</h1>\"
        \"\"\"
        if self.asGroupList:
            raise TypeError(\"cannot use sub() with Regex(asGroupList=True)\")

        if self.asMatch and callable(repl):
            raise TypeError(\"cannot use sub() with a callable with Regex(asMatch=True)\")

        if self.asMatch:

            def pa(tokens):
                return tokens[0].expand(repl)

        else:

            def pa(tokens):
                return self.re.sub(repl, tokens[0])

        return self.add_parse_action(pa)


class QuotedString(Token):
    r\"\"\"
    Token for matching strings that are delimited by quoting characters.

    Defined with the following parameters:

    - ``quote_char`` - string of one or more characters defining the
      quote delimiting string
    - ``esc_char`` - character to re_escape quotes, typically backslash
      (default= ``None``)
    - ``esc_quote`` - special quote sequence to re_escape an embedded quote
      string (such as SQL's ``\"\"`` to re_escape an embedded ``\"``)
      (default= ``None``)
    - ``multiline`` - boolean indicating whether quotes can span
      multiple lines (default= ``False``)
    - ``unquote_results`` - boolean indicating whether the matched text
      should be unquoted (default= ``True``)
    - ``end_quote_char`` - string of one or more characters defining the
      end of the quote delimited string (default= ``None``  => same as
      quote_char)
    - ``convert_whitespace_escapes`` - convert escaped whitespace
      (``'\\t'``, ``'\\n'``, etc.) to actual whitespace
      (default= ``True``)

    Example::

        qs = QuotedString('\"')
        print(qs.search_string('lsjdf \"This is the quote\" sldjf'))
        complex_qs = QuotedString('{{', end_quote_char='}}')
        print(complex_qs.search_string('lsjdf {{This is the \"quote\"}} sldjf'))
        sql_qs = QuotedString('\"', esc_quote='\"\"')
        print(sql_qs.search_string('lsjdf \"This is the quote with \"\"embedded\"\" quotes\" sldjf'))

    prints::

        [['This is the quote']]
        [['This is the \"quote\"']]
        [['This is the quote with \"embedded\" quotes']]
    \"\"\"
    ws_map = ((r\"\\t\", \"\\t\"), (r\"\\n\", \"\\n\"), (r\"\\f\", \"\\f\"), (r\"\\r\", \"\\r\"))

    def __init__(
        self,
        quote_char: str = \"\",
        esc_char: typing.Optional[str] = None,
        esc_quote: typing.Optional[str] = None,
        multiline: bool = False,
        unquote_results: bool = True,
        end_quote_char: typing.Optional[str] = None,
        convert_whitespace_escapes: bool = True,
        *,
        quoteChar: str = \"\",
        escChar: typing.Optional[str] = None,
        escQuote: typing.Optional[str] = None,
        unquoteResults: bool = True,
        endQuoteChar: typing.Optional[str] = None,
        convertWhitespaceEscapes: bool = True,
    ):
        super().__init__()
        escChar = escChar or esc_char
        escQuote = escQuote or esc_quote
        unquoteResults = unquoteResults and unquote_results
        endQuoteChar = endQuoteChar or end_quote_char
        convertWhitespaceEscapes = (
            convertWhitespaceEscapes and convert_whitespace_escapes
        )
        quote_char = quoteChar or quote_char

        # remove white space from quote chars - wont work anyway
        quote_char = quote_char.strip()
        if not quote_char:
            raise ValueError(\"quote_char cannot be the empty string\")

        if endQuoteChar is None:
            endQuoteChar = quote_char
        else:
            endQuoteChar = endQuoteChar.strip()
            if not endQuoteChar:
                raise ValueError(\"endQuoteChar cannot be the empty string\")

        self.quoteChar = quote_char
        self.quoteCharLen = len(quote_char)
        self.firstQuoteChar = quote_char[0]
        self.endQuoteChar = endQuoteChar
        self.endQuoteCharLen = len(endQuoteChar)
        self.escChar = escChar
        self.escQuote = escQuote
        self.unquoteResults = unquoteResults
        self.convertWhitespaceEscapes = convertWhitespaceEscapes

        sep = \"\"
        inner_pattern = \"\"

        if escQuote:
            inner_pattern += r\"{}(?:{})\".format(sep, re.escape(escQuote))
            sep = \"|\"

        if escChar:
            inner_pattern += r\"{}(?:{}.)\".format(sep, re.escape(escChar))
            sep = \"|\"
            self.escCharReplacePattern = re.escape(self.escChar) + \"(.)\"

        if len(self.endQuoteChar) > 1:
            inner_pattern += (
                \"{}(?:\".format(sep)
                + \"|\".join(
                    \"(?:{}(?!{}))\".format(
                        re.escape(self.endQuoteChar[:i]),
                        re.escape(self.endQuoteChar[i:]),
                    )
                    for i in range(len(self.endQuoteChar) - 1, 0, -1)
                )
                + \")\"
            )
            sep = \"|\"

        if multiline:
            self.flags = re.MULTILINE | re.DOTALL
            inner_pattern += r\"{}(?:[^{}{}])\".format(
                sep,
                _escape_regex_range_chars(self.endQuoteChar[0]),
                (_escape_regex_range_chars(escChar) if escChar is not None else \"\"),
            )
        else:
            self.flags = 0
            inner_pattern += r\"{}(?:[^{}\\n\\r{}])\".format(
                sep,
                _escape_regex_range_chars(self.endQuoteChar[0]),
                (_escape_regex_range_chars(escChar) if escChar is not None else \"\"),
            )

        self.pattern = \"\".join(
            [
                re.escape(self.quoteChar),
                \"(?:\",
                inner_pattern,
                \")*\",
                re.escape(self.endQuoteChar),
            ]
        )

        try:
            self.re = re.compile(self.pattern, self.flags)
            self.reString = self.pattern
            self.re_match = self.re.match
        except re.error:
            raise ValueError(
                \"invalid pattern {!r} passed to Regex\".format(self.pattern)
            )

        self.errmsg = \"Expected \" + self.name
        self.mayIndexError = False
        self.mayReturnEmpty = True

    def _generateDefaultName(self):
        if self.quoteChar == self.endQuoteChar and isinstance(self.quoteChar, str_type):
            return \"string enclosed in {!r}\".format(self.quoteChar)

        return \"quoted string, starting with {} ending with {}\".format(
            self.quoteChar, self.endQuoteChar
        )

    def parseImpl(self, instring, loc, doActions=True):
        result = (
            instring[loc] == self.firstQuoteChar
            and self.re_match(instring, loc)
            or None
        )
        if not result:
            raise ParseException(instring, loc, self.errmsg, self)

        loc = result.end()
        ret = result.group()

        if self.unquoteResults:

            # strip off quotes
            ret = ret[self.quoteCharLen : -self.endQuoteCharLen]

            if isinstance(ret, str_type):
                # replace escaped whitespace
                if \"\\\\\" in ret and self.convertWhitespaceEscapes:
                    for wslit, wschar in self.ws_map:
                        ret = ret.replace(wslit, wschar)

                # replace escaped characters
                if self.escChar:
                    ret = re.sub(self.escCharReplacePattern, r\"\\g<1>\", ret)

                # replace escaped quotes
                if self.escQuote:
                    ret = ret.replace(self.escQuote, self.endQuoteChar)

        return loc, ret


class CharsNotIn(Token):
    \"\"\"Token for matching words composed of characters *not* in a given
    set (will include whitespace in matched characters if not listed in
    the provided exclusion set - see example). Defined with string
    containing all disallowed characters, and an optional minimum,
    maximum, and/or exact length.  The default value for ``min`` is
    1 (a minimum value < 1 is not valid); the default values for
    ``max`` and ``exact`` are 0, meaning no maximum or exact
    length restriction.

    Example::

        # define a comma-separated-value as anything that is not a ','
        csv_value = CharsNotIn(',')
        print(delimited_list(csv_value).parse_string(\"dkls,lsdkjf,s12 34,@!#,213\"))

    prints::

        ['dkls', 'lsdkjf', 's12 34', '@!#', '213']
    \"\"\"

    def __init__(
        self,
        not_chars: str = \"\",
        min: int = 1,
        max: int = 0,
        exact: int = 0,
        *,
        notChars: str = \"\",
    ):
        super().__init__()
        self.skipWhitespace = False
        self.notChars = not_chars or notChars
        self.notCharsSet = set(self.notChars)

        if min < 1:
            raise ValueError(
                \"cannot specify a minimum length < 1; use \"
                \"Opt(CharsNotIn()) if zero-length char group is permitted\"
            )

        self.minLen = min

        if max > 0:
            self.maxLen = max
        else:
            self.maxLen = _MAX_INT

        if exact > 0:
            self.maxLen = exact
            self.minLen = exact

        self.errmsg = \"Expected \" + self.name
        self.mayReturnEmpty = self.minLen == 0
        self.mayIndexError = False

    def _generateDefaultName(self):
        not_chars_str = _collapse_string_to_ranges(self.notChars)
        if len(not_chars_str) > 16:
            return \"!W:({}...)\".format(self.notChars[: 16 - 3])
        else:
            return \"!W:({})\".format(self.notChars)

    def parseImpl(self, instring, loc, doActions=True):
        notchars = self.notCharsSet
        if instring[loc] in notchars:
            raise ParseException(instring, loc, self.errmsg, self)

        start = loc
        loc += 1
        maxlen = min(start + self.maxLen, len(instring))
        while loc < maxlen and instring[loc] not in notchars:
            loc += 1

        if loc - start < self.minLen:
            raise ParseException(instring, loc, self.errmsg, self)

        return loc, instring[start:loc]


class White(Token):
    \"\"\"Special matching class for matching whitespace.  Normally,
    whitespace is ignored by pyparsing grammars.  This class is included
    when some whitespace structures are significant.  Define with
    a string containing the whitespace characters to be matched; default
    is ``\" \\\\t\\\\r\\\\n\"``.  Also takes optional ``min``,
    ``max``, and ``exact`` arguments, as defined for the
    :class:`Word` class.
    \"\"\"

    whiteStrs = {
        \" \": \"<SP>\",
        \"\\t\": \"<TAB>\",
        \"\\n\": \"<LF>\",
        \"\\r\": \"<CR>\",
        \"\\f\": \"<FF>\",
        \"\\u00A0\": \"<NBSP>\",
        \"\\u1680\": \"<OGHAM_SPACE_MARK>\",
        \"\\u180E\": \"<MONGOLIAN_VOWEL_SEPARATOR>\",
        \"\\u2000\": \"<EN_QUAD>\",
        \"\\u2001\": \"<EM_QUAD>\",
        \"\\u2002\": \"<EN_SPACE>\",
        \"\\u2003\": \"<EM_SPACE>\",
        \"\\u2004\": \"<THREE-PER-EM_SPACE>\",
        \"\\u2005\": \"<FOUR-PER-EM_SPACE>\",
        \"\\u2006\": \"<SIX-PER-EM_SPACE>\",
        \"\\u2007\": \"<FIGURE_SPACE>\",
        \"\\u2008\": \"<PUNCTUATION_SPACE>\",
        \"\\u2009\": \"<THIN_SPACE>\",
        \"\\u200A\": \"<HAIR_SPACE>\",
        \"\\u200B\": \"<ZERO_WIDTH_SPACE>\",
        \"\\u202F\": \"<NNBSP>\",
        \"\\u205F\": \"<MMSP>\",
        \"\\u3000\": \"<IDEOGRAPHIC_SPACE>\",
    }

    def __init__(self, ws: str = \" \\t\\r\\n\", min: int = 1, max: int = 0, exact: int = 0):
        super().__init__()
        self.matchWhite = ws
        self.set_whitespace_chars(
            \"\".join(c for c in self.whiteStrs if c not in self.matchWhite),
            copy_defaults=True,
        )
        # self.leave_whitespace()
        self.mayReturnEmpty = True
        self.errmsg = \"Expected \" + self.name

        self.minLen = min

        if max > 0:
            self.maxLen = max
        else:
            self.maxLen = _MAX_INT

        if exact > 0:
            self.maxLen = exact
            self.minLen = exact

    def _generateDefaultName(self):
        return \"\".join(White.whiteStrs[c] for c in self.matchWhite)

    def parseImpl(self, instring, loc, doActions=True):
        if instring[loc] not in self.matchWhite:
            raise ParseException(instring, loc, self.errmsg, self)
        start = loc
        loc += 1
        maxloc = start + self.maxLen
        maxloc = min(maxloc, len(instring))
        while loc < maxloc and instring[loc] in self.matchWhite:
            loc += 1

        if loc - start < self.minLen:
            raise ParseException(instring, loc, self.errmsg, self)

        return loc, instring[start:loc]


class PositionToken(Token):
    def __init__(self):
        super().__init__()
        self.mayReturnEmpty = True
        self.mayIndexError = False


class GoToColumn(PositionToken):
    \"\"\"Token to advance to a specific column of input text; useful for
    tabular report scraping.
    \"\"\"

    def __init__(self, colno: int):
        super().__init__()
        self.col = colno

    def preParse(self, instring, loc):
        if col(loc, instring) != self.col:
            instrlen = len(instring)
            if self.ignoreExprs:
                loc = self._skipIgnorables(instring, loc)
            while (
                loc < instrlen
                and instring[loc].isspace()
                and col(loc, instring) != self.col
            ):
                loc += 1
        return loc

    def parseImpl(self, instring, loc, doActions=True):
        thiscol = col(loc, instring)
        if thiscol > self.col:
            raise ParseException(instring, loc, \"Text not in expected column\", self)
        newloc = loc + self.col - thiscol
        ret = instring[loc:newloc]
        return newloc, ret


class LineStart(PositionToken):
    r\"\"\"Matches if current position is at the beginning of a line within
    the parse string

    Example::

        test = '''\\
        AAA this line
        AAA and this line
          AAA but not this one
        B AAA and definitely not this one
        '''

        for t in (LineStart() + 'AAA' + restOfLine).search_string(test):
            print(t)

    prints::

        ['AAA', ' this line']
        ['AAA', ' and this line']

    \"\"\"

    def __init__(self):
        super().__init__()
        self.leave_whitespace()
        self.orig_whiteChars = set() | self.whiteChars
        self.whiteChars.discard(\"\\n\")
        self.skipper = Empty().set_whitespace_chars(self.whiteChars)
        self.errmsg = \"Expected start of line\"

    def preParse(self, instring, loc):
        if loc == 0:
            return loc
        else:
            ret = self.skipper.preParse(instring, loc)
            if \"\\n\" in self.orig_whiteChars:
                while instring[ret : ret + 1] == \"\\n\":
                    ret = self.skipper.preParse(instring, ret + 1)
            return ret

    def parseImpl(self, instring, loc, doActions=True):
        if col(loc, instring) == 1:
            return loc, []
        raise ParseException(instring, loc, self.errmsg, self)


class LineEnd(PositionToken):
    \"\"\"Matches if current position is at the end of a line within the
    parse string
    \"\"\"

    def __init__(self):
        super().__init__()
        self.whiteChars.discard(\"\\n\")
        self.set_whitespace_chars(self.whiteChars, copy_defaults=False)
        self.errmsg = \"Expected end of line\"

    def parseImpl(self, instring, loc, doActions=True):
        if loc < len(instring):
            if instring[loc] == \"\\n\":
                return loc + 1, \"\\n\"
            else:
                raise ParseException(instring, loc, self.errmsg, self)
        elif loc == len(instring):
            return loc + 1, []
        else:
            raise ParseException(instring, loc, self.errmsg, self)


class StringStart(PositionToken):
    \"\"\"Matches if current position is at the beginning of the parse
    string
    \"\"\"

    def __init__(self):
        super().__init__()
        self.errmsg = \"Expected start of text\"

    def parseImpl(self, instring, loc, doActions=True):
        if loc != 0:
            # see if entire string up to here is just whitespace and ignoreables
            if loc != self.preParse(instring, 0):
                raise ParseException(instring, loc, self.errmsg, self)
        return loc, []


class StringEnd(PositionToken):
    \"\"\"
    Matches if current position is at the end of the parse string
    \"\"\"

    def __init__(self):
        super().__init__()
        self.errmsg = \"Expected end of text\"

    def parseImpl(self, instring, loc, doActions=True):
        if loc < len(instring):
            raise ParseException(instring, loc, self.errmsg, self)
        elif loc == len(instring):
            return loc + 1, []
        elif loc > len(instring):
            return loc, []
        else:
            raise ParseException(instring, loc, self.errmsg, self)


class WordStart(PositionToken):
    \"\"\"Matches if the current position is at the beginning of a
    :class:`Word`, and is not preceded by any character in a given
    set of ``word_chars`` (default= ``printables``). To emulate the
    ``\\b`` behavior of regular expressions, use
    ``WordStart(alphanums)``. ``WordStart`` will also match at
    the beginning of the string being parsed, or at the beginning of
    a line.
    \"\"\"

    def __init__(self, word_chars: str = printables, *, wordChars: str = printables):
        wordChars = word_chars if wordChars == printables else wordChars
        super().__init__()
        self.wordChars = set(wordChars)
        self.errmsg = \"Not at the start of a word\"

    def parseImpl(self, instring, loc, doActions=True):
        if loc != 0:
            if (
                instring[loc - 1] in self.wordChars
                or instring[loc] not in self.wordChars
            ):
                raise ParseException(instring, loc, self.errmsg, self)
        return loc, []


class WordEnd(PositionToken):
    \"\"\"Matches if the current position is at the end of a :class:`Word`,
    and is not followed by any character in a given set of ``word_chars``
    (default= ``printables``). To emulate the ``\\b`` behavior of
    regular expressions, use ``WordEnd(alphanums)``. ``WordEnd``
    will also match at the end of the string being parsed, or at the end
    of a line.
    \"\"\"

    def __init__(self, word_chars: str = printables, *, wordChars: str = printables):
        wordChars = word_chars if wordChars == printables else wordChars
        super().__init__()
        self.wordChars = set(wordChars)
        self.skipWhitespace = False
        self.errmsg = \"Not at the end of a word\"

    def parseImpl(self, instring, loc, doActions=True):
        instrlen = len(instring)
        if instrlen > 0 and loc < instrlen:
            if (
                instring[loc] in self.wordChars
                or instring[loc - 1] not in self.wordChars
            ):
                raise ParseException(instring, loc, self.errmsg, self)
        return loc, []


class ParseExpression(ParserElement):
    \"\"\"Abstract subclass of ParserElement, for combining and
    post-processing parsed tokens.
    \"\"\"

    def __init__(self, exprs: typing.Iterable[ParserElement], savelist: bool = False):
        super().__init__(savelist)
        self.exprs: List[ParserElement]
        if isinstance(exprs, _generatorType):
            exprs = list(exprs)

        if isinstance(exprs, str_type):
            self.exprs = [self._literalStringClass(exprs)]
        elif isinstance(exprs, ParserElement):
            self.exprs = [exprs]
        elif isinstance(exprs, Iterable):
            exprs = list(exprs)
            # if sequence of strings provided, wrap with Literal
            if any(isinstance(expr, str_type) for expr in exprs):
                exprs = (
                    self._literalStringClass(e) if isinstance(e, str_type) else e
                    for e in exprs
                )
            self.exprs = list(exprs)
        else:
            try:
                self.exprs = list(exprs)
            except TypeError:
                self.exprs = [exprs]
        self.callPreparse = False

    def recurse(self) -> Sequence[ParserElement]:
        return self.exprs[:]

    def append(self, other) -> ParserElement:
        self.exprs.append(other)
        self._defaultName = None
        return self

    def leave_whitespace(self, recursive: bool = True) -> ParserElement:
        \"\"\"
        Extends ``leave_whitespace`` defined in base class, and also invokes ``leave_whitespace`` on
           all contained expressions.
        \"\"\"
        super().leave_whitespace(recursive)

        if recursive:
            self.exprs = [e.copy() for e in self.exprs]
            for e in self.exprs:
                e.leave_whitespace(recursive)
        return self

    def ignore_whitespace(self, recursive: bool = True) -> ParserElement:
        \"\"\"
        Extends ``ignore_whitespace`` defined in base class, and also invokes ``leave_whitespace`` on
           all contained expressions.
        \"\"\"
        super().ignore_whitespace(recursive)
        if recursive:
            self.exprs = [e.copy() for e in self.exprs]
            for e in self.exprs:
                e.ignore_whitespace(recursive)
        return self

    def ignore(self, other) -> ParserElement:
        if isinstance(other, Suppress):
            if other not in self.ignoreExprs:
                super().ignore(other)
                for e in self.exprs:
                    e.ignore(self.ignoreExprs[-1])
        else:
            super().ignore(other)
            for e in self.exprs:
                e.ignore(self.ignoreExprs[-1])
        return self

    def _generateDefaultName(self):
        return \"{}:({})\".format(self.__class__.__name__, str(self.exprs))

    def streamline(self) -> ParserElement:
        if self.streamlined:
            return self

        super().streamline()

        for e in self.exprs:
            e.streamline()

        # collapse nested :class:`And`'s of the form ``And(And(And(a, b), c), d)`` to ``And(a, b, c, d)``
        # but only if there are no parse actions or resultsNames on the nested And's
        # (likewise for :class:`Or`'s and :class:`MatchFirst`'s)
        if len(self.exprs) == 2:
            other = self.exprs[0]
            if (
                isinstance(other, self.__class__)
                and not other.parseAction
                and other.resultsName is None
                and not other.debug
            ):
                self.exprs = other.exprs[:] + [self.exprs[1]]
                self._defaultName = None
                self.mayReturnEmpty |= other.mayReturnEmpty
                self.mayIndexError |= other.mayIndexError

            other = self.exprs[-1]
            if (
                isinstance(other, self.__class__)
                and not other.parseAction
                and other.resultsName is None
                and not other.debug
            ):
                self.exprs = self.exprs[:-1] + other.exprs[:]
                self._defaultName = None
                self.mayReturnEmpty |= other.mayReturnEmpty
                self.mayIndexError |= other.mayIndexError

        self.errmsg = \"Expected \" + str(self)

        return self

    def validate(self, validateTrace=None) -> None:
        tmp = (validateTrace if validateTrace is not None else [])[:] + [self]
        for e in self.exprs:
            e.validate(tmp)
        self._checkRecursion([])

    def copy(self) -> ParserElement:
        ret = super().copy()
        ret.exprs = [e.copy() for e in self.exprs]
        return ret

    def _setResultsName(self, name, listAllMatches=False):
        if (
            __diag__.warn_ungrouped_named_tokens_in_collection
            and Diagnostics.warn_ungrouped_named_tokens_in_collection
            not in self.suppress_warnings_
        ):
            for e in self.exprs:
                if (
                    isinstance(e, ParserElement)
                    and e.resultsName
                    and Diagnostics.warn_ungrouped_named_tokens_in_collection
                    not in e.suppress_warnings_
                ):
                    warnings.warn(
                        \"{}: setting results name {!r} on {} expression \"
                        \"collides with {!r} on contained expression\".format(
                            \"warn_ungrouped_named_tokens_in_collection\",
                            name,
                            type(self).__name__,
                            e.resultsName,
                        ),
                        stacklevel=3,
                    )

        return super()._setResultsName(name, listAllMatches)

    ignoreWhitespace = ignore_whitespace
    leaveWhitespace = leave_whitespace


class And(ParseExpression):
    \"\"\"
    Requires all given :class:`ParseExpression` s to be found in the given order.
    Expressions may be separated by whitespace.
    May be constructed using the ``'+'`` operator.
    May also be constructed using the ``'-'`` operator, which will
    suppress backtracking.

    Example::

        integer = Word(nums)
        name_expr = Word(alphas)[1, ...]

        expr = And([integer(\"id\"), name_expr(\"name\"), integer(\"age\")])
        # more easily written as:
        expr = integer(\"id\") + name_expr(\"name\") + integer(\"age\")
    \"\"\"

    class _ErrorStop(Empty):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.leave_whitespace()

        def _generateDefaultName(self):
            return \"-\"

    def __init__(
        self, exprs_arg: typing.Iterable[ParserElement], savelist: bool = True
    ):
        exprs: List[ParserElement] = list(exprs_arg)
        if exprs and Ellipsis in exprs:
            tmp = []
            for i, expr in enumerate(exprs):
                if expr is Ellipsis:
                    if i < len(exprs) - 1:
                        skipto_arg: ParserElement = (Empty() + exprs[i + 1]).exprs[-1]
                        tmp.append(SkipTo(skipto_arg)(\"_skipped*\"))
                    else:
                        raise Exception(
                            \"cannot construct And with sequence ending in ...\"
                        )
                else:
                    tmp.append(expr)
            exprs[:] = tmp
        super().__init__(exprs, savelist)
        if self.exprs:
            self.mayReturnEmpty = all(e.mayReturnEmpty for e in self.exprs)
            if not isinstance(self.exprs[0], White):
                self.set_whitespace_chars(
                    self.exprs[0].whiteChars,
                    copy_defaults=self.exprs[0].copyDefaultWhiteChars,
                )
                self.skipWhitespace = self.exprs[0].skipWhitespace
            else:
                self.skipWhitespace = False
        else:
            self.mayReturnEmpty = True
        self.callPreparse = True

    def streamline(self) -> ParserElement:
        # collapse any _PendingSkip's
        if self.exprs:
            if any(
                isinstance(e, ParseExpression)
                and e.exprs
                and isinstance(e.exprs[-1], _PendingSkip)
                for e in self.exprs[:-1]
            ):
                for i, e in enumerate(self.exprs[:-1]):
                    if e is None:
                        continue
                    if (
                        isinstance(e, ParseExpression)
                        and e.exprs
                        and isinstance(e.exprs[-1], _PendingSkip)
                    ):
                        e.exprs[-1] = e.exprs[-1] + self.exprs[i + 1]
                        self.exprs[i + 1] = None
                self.exprs = [e for e in self.exprs if e is not None]

        super().streamline()

        # link any IndentedBlocks to the prior expression
        for prev, cur in zip(self.exprs, self.exprs[1:]):
            # traverse cur or any first embedded expr of cur looking for an IndentedBlock
            # (but watch out for recursive grammar)
            seen = set()
            while cur:
                if id(cur) in seen:
                    break
                seen.add(id(cur))
                if isinstance(cur, IndentedBlock):
                    prev.add_parse_action(
                        lambda s, l, t, cur_=cur: setattr(
                            cur_, \"parent_anchor\", col(l, s)
                        )
                    )
                    break
                subs = cur.recurse()
                cur = next(iter(subs), None)

        self.mayReturnEmpty = all(e.mayReturnEmpty for e in self.exprs)
        return self

    def parseImpl(self, instring, loc, doActions=True):
        # pass False as callPreParse arg to _parse for first element, since we already
        # pre-parsed the string as part of our And pre-parsing
        loc, resultlist = self.exprs[0]._parse(
            instring, loc, doActions, callPreParse=False
        )
        errorStop = False
        for e in self.exprs[1:]:
            # if isinstance(e, And._ErrorStop):
            if type(e) is And._ErrorStop:
                errorStop = True
                continue
            if errorStop:
                try:
                    loc, exprtokens = e._parse(instring, loc, doActions)
                except ParseSyntaxException:
                    raise
                except ParseBaseException as pe:
                    pe.__traceback__ = None
                    raise ParseSyntaxException._from_exception(pe)
                except IndexError:
                    raise ParseSyntaxException(
                        instring, len(instring), self.errmsg, self
                    )
            else:
                loc, exprtokens = e._parse(instring, loc, doActions)
            if exprtokens or exprtokens.haskeys():
                resultlist += exprtokens
        return loc, resultlist

    def __iadd__(self, other):
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        return self.append(other)  # And([self, other])

    def _checkRecursion(self, parseElementList):
        subRecCheckList = parseElementList[:] + [self]
        for e in self.exprs:
            e._checkRecursion(subRecCheckList)
            if not e.mayReturnEmpty:
                break

    def _generateDefaultName(self):
        inner = \" \".join(str(e) for e in self.exprs)
        # strip off redundant inner {}'s
        while len(inner) > 1 and inner[0 :: len(inner) - 1] == \"{}\":
            inner = inner[1:-1]
        return \"{\" + inner + \"}\"


class Or(ParseExpression):
    \"\"\"Requires that at least one :class:`ParseExpression` is found. If
    two expressions match, the expression that matches the longest
    string will be used. May be constructed using the ``'^'``
    operator.

    Example::

        # construct Or using '^' operator

        number = Word(nums) ^ Combine(Word(nums) + '.' + Word(nums))
        print(number.search_string(\"123 3.1416 789\"))

    prints::

        [['123'], ['3.1416'], ['789']]
    \"\"\"

    def __init__(self, exprs: typing.Iterable[ParserElement], savelist: bool = False):
        super().__init__(exprs, savelist)
        if self.exprs:
            self.mayReturnEmpty = any(e.mayReturnEmpty for e in self.exprs)
            self.skipWhitespace = all(e.skipWhitespace for e in self.exprs)
        else:
            self.mayReturnEmpty = True

    def streamline(self) -> ParserElement:
        super().streamline()
        if self.exprs:
            self.mayReturnEmpty = any(e.mayReturnEmpty for e in self.exprs)
            self.saveAsList = any(e.saveAsList for e in self.exprs)
            self.skipWhitespace = all(
                e.skipWhitespace and not isinstance(e, White) for e in self.exprs
            )
        else:
            self.saveAsList = False
        return self

    def parseImpl(self, instring, loc, doActions=True):
        maxExcLoc = -1
        maxException = None
        matches = []
        fatals = []
        if all(e.callPreparse for e in self.exprs):
            loc = self.preParse(instring, loc)
        for e in self.exprs:
            try:
                loc2 = e.try_parse(instring, loc, raise_fatal=True)
            except ParseFatalException as pfe:
                pfe.__traceback__ = None
                pfe.parserElement = e
                fatals.append(pfe)
                maxException = None
                maxExcLoc = -1
            except ParseException as err:
                if not fatals:
                    err.__traceback__ = None
                    if err.loc > maxExcLoc:
                        maxException = err
                        maxExcLoc = err.loc
            except IndexError:
                if len(instring) > maxExcLoc:
                    maxException = ParseException(
                        instring, len(instring), e.errmsg, self
                    )
                    maxExcLoc = len(instring)
            else:
                # save match among all matches, to retry longest to shortest
                matches.append((loc2, e))

        if matches:
            # re-evaluate all matches in descending order of length of match, in case attached actions
            # might change whether or how much they match of the input.
            matches.sort(key=itemgetter(0), reverse=True)

            if not doActions:
                # no further conditions or parse actions to change the selection of
                # alternative, so the first match will be the best match
                best_expr = matches[0][1]
                return best_expr._parse(instring, loc, doActions)

            longest = -1, None
            for loc1, expr1 in matches:
                if loc1 <= longest[0]:
                    # already have a longer match than this one will deliver, we are done
                    return longest

                try:
                    loc2, toks = expr1._parse(instring, loc, doActions)
                except ParseException as err:
                    err.__traceback__ = None
                    if err.loc > maxExcLoc:
                        maxException = err
                        maxExcLoc = err.loc
                else:
                    if loc2 >= loc1:
                        return loc2, toks
                    # didn't match as much as before
                    elif loc2 > longest[0]:
                        longest = loc2, toks

            if longest != (-1, None):
                return longest

        if fatals:
            if len(fatals) > 1:
                fatals.sort(key=lambda e: -e.loc)
                if fatals[0].loc == fatals[1].loc:
                    fatals.sort(key=lambda e: (-e.loc, -len(str(e.parserElement))))
            max_fatal = fatals[0]
            raise max_fatal

        if maxException is not None:
            maxException.msg = self.errmsg
            raise maxException
        else:
            raise ParseException(
                instring, loc, \"no defined alternatives to match\", self
            )

    def __ixor__(self, other):
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        return self.append(other)  # Or([self, other])

    def _generateDefaultName(self):
        return \"{\" + \" ^ \".join(str(e) for e in self.exprs) + \"}\"

    def _setResultsName(self, name, listAllMatches=False):
        if (
            __diag__.warn_multiple_tokens_in_named_alternation
            and Diagnostics.warn_multiple_tokens_in_named_alternation
            not in self.suppress_warnings_
        ):
            if any(
                isinstance(e, And)
                and Diagnostics.warn_multiple_tokens_in_named_alternation
                not in e.suppress_warnings_
                for e in self.exprs
            ):
                warnings.warn(
                    \"{}: setting results name {!r} on {} expression \"
                    \"will return a list of all parsed tokens in an And alternative, \"
                    \"in prior versions only the first token was returned; enclose \"
                    \"contained argument in Group\".format(
                        \"warn_multiple_tokens_in_named_alternation\",
                        name,
                        type(self).__name__,
                    ),
                    stacklevel=3,
                )

        return super()._setResultsName(name, listAllMatches)


class MatchFirst(ParseExpression):
    \"\"\"Requires that at least one :class:`ParseExpression` is found. If
    more than one expression matches, the first one listed is the one that will
    match. May be constructed using the ``'|'`` operator.

    Example::

        # construct MatchFirst using '|' operator

        # watch the order of expressions to match
        number = Word(nums) | Combine(Word(nums) + '.' + Word(nums))
        print(number.search_string(\"123 3.1416 789\")) #  Fail! -> [['123'], ['3'], ['1416'], ['789']]

        # put more selective expression first
        number = Combine(Word(nums) + '.' + Word(nums)) | Word(nums)
        print(number.search_string(\"123 3.1416 789\")) #  Better -> [['123'], ['3.1416'], ['789']]
    \"\"\"

    def __init__(self, exprs: typing.Iterable[ParserElement], savelist: bool = False):
        super().__init__(exprs, savelist)
        if self.exprs:
            self.mayReturnEmpty = any(e.mayReturnEmpty for e in self.exprs)
            self.skipWhitespace = all(e.skipWhitespace for e in self.exprs)
        else:
            self.mayReturnEmpty = True

    def streamline(self) -> ParserElement:
        if self.streamlined:
            return self

        super().streamline()
        if self.exprs:
            self.saveAsList = any(e.saveAsList for e in self.exprs)
            self.mayReturnEmpty = any(e.mayReturnEmpty for e in self.exprs)
            self.skipWhitespace = all(
                e.skipWhitespace and not isinstance(e, White) for e in self.exprs
            )
        else:
            self.saveAsList = False
            self.mayReturnEmpty = True
        return self

    def parseImpl(self, instring, loc, doActions=True):
        maxExcLoc = -1
        maxException = None

        for e in self.exprs:
            try:
                return e._parse(
                    instring,
                    loc,
                    doActions,
                )
            except ParseFatalException as pfe:
                pfe.__traceback__ = None
                pfe.parserElement = e
                raise
            except ParseException as err:
                if err.loc > maxExcLoc:
                    maxException = err
                    maxExcLoc = err.loc
            except IndexError:
                if len(instring) > maxExcLoc:
                    maxException = ParseException(
                        instring, len(instring), e.errmsg, self
                    )
                    maxExcLoc = len(instring)

        if maxException is not None:
            maxException.msg = self.errmsg
            raise maxException
        else:
            raise ParseException(
                instring, loc, \"no defined alternatives to match\", self
            )

    def __ior__(self, other):
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        return self.append(other)  # MatchFirst([self, other])

    def _generateDefaultName(self):
        return \"{\" + \" | \".join(str(e) for e in self.exprs) + \"}\"

    def _setResultsName(self, name, listAllMatches=False):
        if (
            __diag__.warn_multiple_tokens_in_named_alternation
            and Diagnostics.warn_multiple_tokens_in_named_alternation
            not in self.suppress_warnings_
        ):
            if any(
                isinstance(e, And)
                and Diagnostics.warn_multiple_tokens_in_named_alternation
                not in e.suppress_warnings_
                for e in self.exprs
            ):
                warnings.warn(
                    \"{}: setting results name {!r} on {} expression \"
                    \"will return a list of all parsed tokens in an And alternative, \"
                    \"in prior versions only the first token was returned; enclose \"
                    \"contained argument in Group\".format(
                        \"warn_multiple_tokens_in_named_alternation\",
                        name,
                        type(self).__name__,
                    ),
                    stacklevel=3,
                )

        return super()._setResultsName(name, listAllMatches)


class Each(ParseExpression):
    \"\"\"Requires all given :class:`ParseExpression` s to be found, but in
    any order. Expressions may be separated by whitespace.

    May be constructed using the ``'&'`` operator.

    Example::

        color = one_of(\"RED ORANGE YELLOW GREEN BLUE PURPLE BLACK WHITE BROWN\")
        shape_type = one_of(\"SQUARE CIRCLE TRIANGLE STAR HEXAGON OCTAGON\")
        integer = Word(nums)
        shape_attr = \"shape:\" + shape_type(\"shape\")
        posn_attr = \"posn:\" + Group(integer(\"x\") + ',' + integer(\"y\"))(\"posn\")
        color_attr = \"color:\" + color(\"color\")
        size_attr = \"size:\" + integer(\"size\")

        # use Each (using operator '&') to accept attributes in any order
        # (shape and posn are required, color and size are optional)
        shape_spec = shape_attr & posn_attr & Opt(color_attr) & Opt(size_attr)

        shape_spec.run_tests('''
            shape: SQUARE color: BLACK posn: 100, 120
            shape: CIRCLE size: 50 color: BLUE posn: 50,80
            color:GREEN size:20 shape:TRIANGLE posn:20,40
            '''
            )

    prints::

        shape: SQUARE color: BLACK posn: 100, 120
        ['shape:', 'SQUARE', 'color:', 'BLACK', 'posn:', ['100', ',', '120']]
        - color: BLACK
        - posn: ['100', ',', '120']
          - x: 100
          - y: 120
        - shape: SQUARE


        shape: CIRCLE size: 50 color: BLUE posn: 50,80
        ['shape:', 'CIRCLE', 'size:', '50', 'color:', 'BLUE', 'posn:', ['50', ',', '80']]
        - color: BLUE
        - posn: ['50', ',', '80']
          - x: 50
          - y: 80
        - shape: CIRCLE
        - size: 50


        color: GREEN size: 20 shape: TRIANGLE posn: 20,40
        ['color:', 'GREEN', 'size:', '20', 'shape:', 'TRIANGLE', 'posn:', ['20', ',', '40']]
        - color: GREEN
        - posn: ['20', ',', '40']
          - x: 20
          - y: 40
        - shape: TRIANGLE
        - size: 20
    \"\"\"

    def __init__(self, exprs: typing.Iterable[ParserElement], savelist: bool = True):
        super().__init__(exprs, savelist)
        if self.exprs:
            self.mayReturnEmpty = all(e.mayReturnEmpty for e in self.exprs)
        else:
            self.mayReturnEmpty = True
        self.skipWhitespace = True
        self.initExprGroups = True
        self.saveAsList = True

    def streamline(self) -> ParserElement:
        super().streamline()
        if self.exprs:
            self.mayReturnEmpty = all(e.mayReturnEmpty for e in self.exprs)
        else:
            self.mayReturnEmpty = True
        return self

    def parseImpl(self, instring, loc, doActions=True):
        if self.initExprGroups:
            self.opt1map = dict(
                (id(e.expr), e) for e in self.exprs if isinstance(e, Opt)
            )
            opt1 = [e.expr for e in self.exprs if isinstance(e, Opt)]
            opt2 = [
                e
                for e in self.exprs
                if e.mayReturnEmpty and not isinstance(e, (Opt, Regex, ZeroOrMore))
            ]
            self.optionals = opt1 + opt2
            self.multioptionals = [
                e.expr.set_results_name(e.resultsName, list_all_matches=True)
                for e in self.exprs
                if isinstance(e, _MultipleMatch)
            ]
            self.multirequired = [
                e.expr.set_results_name(e.resultsName, list_all_matches=True)
                for e in self.exprs
                if isinstance(e, OneOrMore)
            ]
            self.required = [
                e for e in self.exprs if not isinstance(e, (Opt, ZeroOrMore, OneOrMore))
            ]
            self.required += self.multirequired
            self.initExprGroups = False

        tmpLoc = loc
        tmpReqd = self.required[:]
        tmpOpt = self.optionals[:]
        multis = self.multioptionals[:]
        matchOrder = []

        keepMatching = True
        failed = []
        fatals = []
        while keepMatching:
            tmpExprs = tmpReqd + tmpOpt + multis
            failed.clear()
            fatals.clear()
            for e in tmpExprs:
                try:
                    tmpLoc = e.try_parse(instring, tmpLoc, raise_fatal=True)
                except ParseFatalException as pfe:
                    pfe.__traceback__ = None
                    pfe.parserElement = e
                    fatals.append(pfe)
                    failed.append(e)
                except ParseException:
                    failed.append(e)
                else:
                    matchOrder.append(self.opt1map.get(id(e), e))
                    if e in tmpReqd:
                        tmpReqd.remove(e)
                    elif e in tmpOpt:
                        tmpOpt.remove(e)
            if len(failed) == len(tmpExprs):
                keepMatching = False

        # look for any ParseFatalExceptions
        if fatals:
            if len(fatals) > 1:
                fatals.sort(key=lambda e: -e.loc)
                if fatals[0].loc == fatals[1].loc:
                    fatals.sort(key=lambda e: (-e.loc, -len(str(e.parserElement))))
            max_fatal = fatals[0]
            raise max_fatal

        if tmpReqd:
            missing = \", \".join([str(e) for e in tmpReqd])
            raise ParseException(
                instring,
                loc,
                \"Missing one or more required elements ({})\".format(missing),
            )

        # add any unmatched Opts, in case they have default values defined
        matchOrder += [e for e in self.exprs if isinstance(e, Opt) and e.expr in tmpOpt]

        total_results = ParseResults([])
        for e in matchOrder:
            loc, results = e._parse(instring, loc, doActions)
            total_results += results

        return loc, total_results

    def _generateDefaultName(self):
        return \"{\" + \" & \".join(str(e) for e in self.exprs) + \"}\"


class ParseElementEnhance(ParserElement):
    \"\"\"Abstract subclass of :class:`ParserElement`, for combining and
    post-processing parsed tokens.
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str], savelist: bool = False):
        super().__init__(savelist)
        if isinstance(expr, str_type):
            if issubclass(self._literalStringClass, Token):
                expr = self._literalStringClass(expr)
            elif issubclass(type(self), self._literalStringClass):
                expr = Literal(expr)
            else:
                expr = self._literalStringClass(Literal(expr))
        self.expr = expr
        if expr is not None:
            self.mayIndexError = expr.mayIndexError
            self.mayReturnEmpty = expr.mayReturnEmpty
            self.set_whitespace_chars(
                expr.whiteChars, copy_defaults=expr.copyDefaultWhiteChars
            )
            self.skipWhitespace = expr.skipWhitespace
            self.saveAsList = expr.saveAsList
            self.callPreparse = expr.callPreparse
            self.ignoreExprs.extend(expr.ignoreExprs)

    def recurse(self) -> Sequence[ParserElement]:
        return [self.expr] if self.expr is not None else []

    def parseImpl(self, instring, loc, doActions=True):
        if self.expr is not None:
            return self.expr._parse(instring, loc, doActions, callPreParse=False)
        else:
            raise ParseException(instring, loc, \"No expression defined\", self)

    def leave_whitespace(self, recursive: bool = True) -> ParserElement:
        super().leave_whitespace(recursive)

        if recursive:
            self.expr = self.expr.copy()
            if self.expr is not None:
                self.expr.leave_whitespace(recursive)
        return self

    def ignore_whitespace(self, recursive: bool = True) -> ParserElement:
        super().ignore_whitespace(recursive)

        if recursive:
            self.expr = self.expr.copy()
            if self.expr is not None:
                self.expr.ignore_whitespace(recursive)
        return self

    def ignore(self, other) -> ParserElement:
        if isinstance(other, Suppress):
            if other not in self.ignoreExprs:
                super().ignore(other)
                if self.expr is not None:
                    self.expr.ignore(self.ignoreExprs[-1])
        else:
            super().ignore(other)
            if self.expr is not None:
                self.expr.ignore(self.ignoreExprs[-1])
        return self

    def streamline(self) -> ParserElement:
        super().streamline()
        if self.expr is not None:
            self.expr.streamline()
        return self

    def _checkRecursion(self, parseElementList):
        if self in parseElementList:
            raise RecursiveGrammarException(parseElementList + [self])
        subRecCheckList = parseElementList[:] + [self]
        if self.expr is not None:
            self.expr._checkRecursion(subRecCheckList)

    def validate(self, validateTrace=None) -> None:
        if validateTrace is None:
            validateTrace = []
        tmp = validateTrace[:] + [self]
        if self.expr is not None:
            self.expr.validate(tmp)
        self._checkRecursion([])

    def _generateDefaultName(self):
        return \"{}:({})\".format(self.__class__.__name__, str(self.expr))

    ignoreWhitespace = ignore_whitespace
    leaveWhitespace = leave_whitespace


class IndentedBlock(ParseElementEnhance):
    \"\"\"
    Expression to match one or more expressions at a given indentation level.
    Useful for parsing text where structure is implied by indentation (like Python source code).
    \"\"\"

    class _Indent(Empty):
        def __init__(self, ref_col: int):
            super().__init__()
            self.errmsg = \"expected indent at column {}\".format(ref_col)
            self.add_condition(lambda s, l, t: col(l, s) == ref_col)

    class _IndentGreater(Empty):
        def __init__(self, ref_col: int):
            super().__init__()
            self.errmsg = \"expected indent at column greater than {}\".format(ref_col)
            self.add_condition(lambda s, l, t: col(l, s) > ref_col)

    def __init__(
        self, expr: ParserElement, *, recursive: bool = False, grouped: bool = True
    ):
        super().__init__(expr, savelist=True)
        # if recursive:
        #     raise NotImplementedError(\"IndentedBlock with recursive is not implemented\")
        self._recursive = recursive
        self._grouped = grouped
        self.parent_anchor = 1

    def parseImpl(self, instring, loc, doActions=True):
        # advance parse position to non-whitespace by using an Empty()
        # this should be the column to be used for all subsequent indented lines
        anchor_loc = Empty().preParse(instring, loc)

        # see if self.expr matches at the current location - if not it will raise an exception
        # and no further work is necessary
        self.expr.try_parse(instring, anchor_loc, doActions)

        indent_col = col(anchor_loc, instring)
        peer_detect_expr = self._Indent(indent_col)

        inner_expr = Empty() + peer_detect_expr + self.expr
        if self._recursive:
            sub_indent = self._IndentGreater(indent_col)
            nested_block = IndentedBlock(
                self.expr, recursive=self._recursive, grouped=self._grouped
            )
            nested_block.set_debug(self.debug)
            nested_block.parent_anchor = indent_col
            inner_expr += Opt(sub_indent + nested_block)

        inner_expr.set_name(f\"inner {hex(id(inner_expr))[-4:].upper()}@{indent_col}\")
        block = OneOrMore(inner_expr)

        trailing_undent = self._Indent(self.parent_anchor) | StringEnd()

        if self._grouped:
            wrapper = Group
        else:
            wrapper = lambda expr: expr
        return (wrapper(block) + Optional(trailing_undent)).parseImpl(
            instring, anchor_loc, doActions
        )


class AtStringStart(ParseElementEnhance):
    \"\"\"Matches if expression matches at the beginning of the parse
    string::

        AtStringStart(Word(nums)).parse_string(\"123\")
        # prints [\"123\"]

        AtStringStart(Word(nums)).parse_string(\"    123\")
        # raises ParseException
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str]):
        super().__init__(expr)
        self.callPreparse = False

    def parseImpl(self, instring, loc, doActions=True):
        if loc != 0:
            raise ParseException(instring, loc, \"not found at string start\")
        return super().parseImpl(instring, loc, doActions)


class AtLineStart(ParseElementEnhance):
    r\"\"\"Matches if an expression matches at the beginning of a line within
    the parse string

    Example::

        test = '''\\
        AAA this line
        AAA and this line
          AAA but not this one
        B AAA and definitely not this one
        '''

        for t in (AtLineStart('AAA') + restOfLine).search_string(test):
            print(t)

    prints::

        ['AAA', ' this line']
        ['AAA', ' and this line']

    \"\"\"

    def __init__(self, expr: Union[ParserElement, str]):
        super().__init__(expr)
        self.callPreparse = False

    def parseImpl(self, instring, loc, doActions=True):
        if col(loc, instring) != 1:
            raise ParseException(instring, loc, \"not found at line start\")
        return super().parseImpl(instring, loc, doActions)


class FollowedBy(ParseElementEnhance):
    \"\"\"Lookahead matching of the given parse expression.
    ``FollowedBy`` does *not* advance the parsing position within
    the input string, it only verifies that the specified parse
    expression matches at the current position.  ``FollowedBy``
    always returns a null token list. If any results names are defined
    in the lookahead expression, those *will* be returned for access by
    name.

    Example::

        # use FollowedBy to match a label only if it is followed by a ':'
        data_word = Word(alphas)
        label = data_word + FollowedBy(':')
        attr_expr = Group(label + Suppress(':') + OneOrMore(data_word, stop_on=label).set_parse_action(' '.join))

        attr_expr[1, ...].parse_string(\"shape: SQUARE color: BLACK posn: upper left\").pprint()

    prints::

        [['shape', 'SQUARE'], ['color', 'BLACK'], ['posn', 'upper left']]
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str]):
        super().__init__(expr)
        self.mayReturnEmpty = True

    def parseImpl(self, instring, loc, doActions=True):
        # by using self._expr.parse and deleting the contents of the returned ParseResults list
        # we keep any named results that were defined in the FollowedBy expression
        _, ret = self.expr._parse(instring, loc, doActions=doActions)
        del ret[:]

        return loc, ret


class PrecededBy(ParseElementEnhance):
    \"\"\"Lookbehind matching of the given parse expression.
    ``PrecededBy`` does not advance the parsing position within the
    input string, it only verifies that the specified parse expression
    matches prior to the current position.  ``PrecededBy`` always
    returns a null token list, but if a results name is defined on the
    given expression, it is returned.

    Parameters:

    - expr - expression that must match prior to the current parse
      location
    - retreat - (default= ``None``) - (int) maximum number of characters
      to lookbehind prior to the current parse location

    If the lookbehind expression is a string, :class:`Literal`,
    :class:`Keyword`, or a :class:`Word` or :class:`CharsNotIn`
    with a specified exact or maximum length, then the retreat
    parameter is not required. Otherwise, retreat must be specified to
    give a maximum number of characters to look back from
    the current parse position for a lookbehind match.

    Example::

        # VB-style variable names with type prefixes
        int_var = PrecededBy(\"#\") + pyparsing_common.identifier
        str_var = PrecededBy(\"$\") + pyparsing_common.identifier

    \"\"\"

    def __init__(
        self, expr: Union[ParserElement, str], retreat: typing.Optional[int] = None
    ):
        super().__init__(expr)
        self.expr = self.expr().leave_whitespace()
        self.mayReturnEmpty = True
        self.mayIndexError = False
        self.exact = False
        if isinstance(expr, str_type):
            retreat = len(expr)
            self.exact = True
        elif isinstance(expr, (Literal, Keyword)):
            retreat = expr.matchLen
            self.exact = True
        elif isinstance(expr, (Word, CharsNotIn)) and expr.maxLen != _MAX_INT:
            retreat = expr.maxLen
            self.exact = True
        elif isinstance(expr, PositionToken):
            retreat = 0
            self.exact = True
        self.retreat = retreat
        self.errmsg = \"not preceded by \" + str(expr)
        self.skipWhitespace = False
        self.parseAction.append(lambda s, l, t: t.__delitem__(slice(None, None)))

    def parseImpl(self, instring, loc=0, doActions=True):
        if self.exact:
            if loc < self.retreat:
                raise ParseException(instring, loc, self.errmsg)
            start = loc - self.retreat
            _, ret = self.expr._parse(instring, start)
        else:
            # retreat specified a maximum lookbehind window, iterate
            test_expr = self.expr + StringEnd()
            instring_slice = instring[max(0, loc - self.retreat) : loc]
            last_expr = ParseException(instring, loc, self.errmsg)
            for offset in range(1, min(loc, self.retreat + 1) + 1):
                try:
                    # print('trying', offset, instring_slice, repr(instring_slice[loc - offset:]))
                    _, ret = test_expr._parse(
                        instring_slice, len(instring_slice) - offset
                    )
                except ParseBaseException as pbe:
                    last_expr = pbe
                else:
                    break
            else:
                raise last_expr
        return loc, ret


class Located(ParseElementEnhance):
    \"\"\"
    Decorates a returned token with its starting and ending
    locations in the input string.

    This helper adds the following results names:

    - ``locn_start`` - location where matched expression begins
    - ``locn_end`` - location where matched expression ends
    - ``value`` - the actual parsed results

    Be careful if the input text contains ``<TAB>`` characters, you
    may want to call :class:`ParserElement.parse_with_tabs`

    Example::

        wd = Word(alphas)
        for match in Located(wd).search_string(\"ljsdf123lksdjjf123lkkjj1222\"):
            print(match)

    prints::

        [0, ['ljsdf'], 5]
        [8, ['lksdjjf'], 15]
        [18, ['lkkjj'], 23]

    \"\"\"

    def parseImpl(self, instring, loc, doActions=True):
        start = loc
        loc, tokens = self.expr._parse(instring, start, doActions, callPreParse=False)
        ret_tokens = ParseResults([start, tokens, loc])
        ret_tokens[\"locn_start\"] = start
        ret_tokens[\"value\"] = tokens
        ret_tokens[\"locn_end\"] = loc
        if self.resultsName:
            # must return as a list, so that the name will be attached to the complete group
            return loc, [ret_tokens]
        else:
            return loc, ret_tokens


class NotAny(ParseElementEnhance):
    \"\"\"
    Lookahead to disallow matching with the given parse expression.
    ``NotAny`` does *not* advance the parsing position within the
    input string, it only verifies that the specified parse expression
    does *not* match at the current position.  Also, ``NotAny`` does
    *not* skip over leading whitespace. ``NotAny`` always returns
    a null token list.  May be constructed using the ``'~'`` operator.

    Example::

        AND, OR, NOT = map(CaselessKeyword, \"AND OR NOT\".split())

        # take care not to mistake keywords for identifiers
        ident = ~(AND | OR | NOT) + Word(alphas)
        boolean_term = Opt(NOT) + ident

        # very crude boolean expression - to support parenthesis groups and
        # operation hierarchy, use infix_notation
        boolean_expr = boolean_term + ((AND | OR) + boolean_term)[...]

        # integers that are followed by \".\" are actually floats
        integer = Word(nums) + ~Char(\".\")
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str]):
        super().__init__(expr)
        # do NOT use self.leave_whitespace(), don't want to propagate to exprs
        # self.leave_whitespace()
        self.skipWhitespace = False

        self.mayReturnEmpty = True
        self.errmsg = \"Found unwanted token, \" + str(self.expr)

    def parseImpl(self, instring, loc, doActions=True):
        if self.expr.can_parse_next(instring, loc):
            raise ParseException(instring, loc, self.errmsg, self)
        return loc, []

    def _generateDefaultName(self):
        return \"~{\" + str(self.expr) + \"}\"


class _MultipleMatch(ParseElementEnhance):
    def __init__(
        self,
        expr: ParserElement,
        stop_on: typing.Optional[Union[ParserElement, str]] = None,
        *,
        stopOn: typing.Optional[Union[ParserElement, str]] = None,
    ):
        super().__init__(expr)
        stopOn = stopOn or stop_on
        self.saveAsList = True
        ender = stopOn
        if isinstance(ender, str_type):
            ender = self._literalStringClass(ender)
        self.stopOn(ender)

    def stopOn(self, ender) -> ParserElement:
        if isinstance(ender, str_type):
            ender = self._literalStringClass(ender)
        self.not_ender = ~ender if ender is not None else None
        return self

    def parseImpl(self, instring, loc, doActions=True):
        self_expr_parse = self.expr._parse
        self_skip_ignorables = self._skipIgnorables
        check_ender = self.not_ender is not None
        if check_ender:
            try_not_ender = self.not_ender.tryParse

        # must be at least one (but first see if we are the stopOn sentinel;
        # if so, fail)
        if check_ender:
            try_not_ender(instring, loc)
        loc, tokens = self_expr_parse(instring, loc, doActions)
        try:
            hasIgnoreExprs = not not self.ignoreExprs
            while 1:
                if check_ender:
                    try_not_ender(instring, loc)
                if hasIgnoreExprs:
                    preloc = self_skip_ignorables(instring, loc)
                else:
                    preloc = loc
                loc, tmptokens = self_expr_parse(instring, preloc, doActions)
                if tmptokens or tmptokens.haskeys():
                    tokens += tmptokens
        except (ParseException, IndexError):
            pass

        return loc, tokens

    def _setResultsName(self, name, listAllMatches=False):
        if (
            __diag__.warn_ungrouped_named_tokens_in_collection
            and Diagnostics.warn_ungrouped_named_tokens_in_collection
            not in self.suppress_warnings_
        ):
            for e in [self.expr] + self.expr.recurse():
                if (
                    isinstance(e, ParserElement)
                    and e.resultsName
                    and Diagnostics.warn_ungrouped_named_tokens_in_collection
                    not in e.suppress_warnings_
                ):
                    warnings.warn(
                        \"{}: setting results name {!r} on {} expression \"
                        \"collides with {!r} on contained expression\".format(
                            \"warn_ungrouped_named_tokens_in_collection\",
                            name,
                            type(self).__name__,
                            e.resultsName,
                        ),
                        stacklevel=3,
                    )

        return super()._setResultsName(name, listAllMatches)


class OneOrMore(_MultipleMatch):
    \"\"\"
    Repetition of one or more of the given expression.

    Parameters:
    - expr - expression that must match one or more times
    - stop_on - (default= ``None``) - expression for a terminating sentinel
         (only required if the sentinel would ordinarily match the repetition
         expression)

    Example::

        data_word = Word(alphas)
        label = data_word + FollowedBy(':')
        attr_expr = Group(label + Suppress(':') + OneOrMore(data_word).set_parse_action(' '.join))

        text = \"shape: SQUARE posn: upper left color: BLACK\"
        attr_expr[1, ...].parse_string(text).pprint()  # Fail! read 'color' as data instead of next label -> [['shape', 'SQUARE color']]

        # use stop_on attribute for OneOrMore to avoid reading label string as part of the data
        attr_expr = Group(label + Suppress(':') + OneOrMore(data_word, stop_on=label).set_parse_action(' '.join))
        OneOrMore(attr_expr).parse_string(text).pprint() # Better -> [['shape', 'SQUARE'], ['posn', 'upper left'], ['color', 'BLACK']]

        # could also be written as
        (attr_expr * (1,)).parse_string(text).pprint()
    \"\"\"

    def _generateDefaultName(self):
        return \"{\" + str(self.expr) + \"}...\"


class ZeroOrMore(_MultipleMatch):
    \"\"\"
    Optional repetition of zero or more of the given expression.

    Parameters:
    - ``expr`` - expression that must match zero or more times
    - ``stop_on`` - expression for a terminating sentinel
      (only required if the sentinel would ordinarily match the repetition
      expression) - (default= ``None``)

    Example: similar to :class:`OneOrMore`
    \"\"\"

    def __init__(
        self,
        expr: ParserElement,
        stop_on: typing.Optional[Union[ParserElement, str]] = None,
        *,
        stopOn: typing.Optional[Union[ParserElement, str]] = None,
    ):
        super().__init__(expr, stopOn=stopOn or stop_on)
        self.mayReturnEmpty = True

    def parseImpl(self, instring, loc, doActions=True):
        try:
            return super().parseImpl(instring, loc, doActions)
        except (ParseException, IndexError):
            return loc, ParseResults([], name=self.resultsName)

    def _generateDefaultName(self):
        return \"[\" + str(self.expr) + \"]...\"


class _NullToken:
    def __bool__(self):
        return False

    def __str__(self):
        return \"\"


class Opt(ParseElementEnhance):
    \"\"\"
    Optional matching of the given expression.

    Parameters:
    - ``expr`` - expression that must match zero or more times
    - ``default`` (optional) - value to be returned if the optional expression is not found.

    Example::

        # US postal code can be a 5-digit zip, plus optional 4-digit qualifier
        zip = Combine(Word(nums, exact=5) + Opt('-' + Word(nums, exact=4)))
        zip.run_tests('''
            # traditional ZIP code
            12345

            # ZIP+4 form
            12101-0001

            # invalid ZIP
            98765-
            ''')

    prints::

        # traditional ZIP code
        12345
        ['12345']

        # ZIP+4 form
        12101-0001
        ['12101-0001']

        # invalid ZIP
        98765-
             ^
        FAIL: Expected end of text (at char 5), (line:1, col:6)
    \"\"\"

    __optionalNotMatched = _NullToken()

    def __init__(
        self, expr: Union[ParserElement, str], default: Any = __optionalNotMatched
    ):
        super().__init__(expr, savelist=False)
        self.saveAsList = self.expr.saveAsList
        self.defaultValue = default
        self.mayReturnEmpty = True

    def parseImpl(self, instring, loc, doActions=True):
        self_expr = self.expr
        try:
            loc, tokens = self_expr._parse(instring, loc, doActions, callPreParse=False)
        except (ParseException, IndexError):
            default_value = self.defaultValue
            if default_value is not self.__optionalNotMatched:
                if self_expr.resultsName:
                    tokens = ParseResults([default_value])
                    tokens[self_expr.resultsName] = default_value
                else:
                    tokens = [default_value]
            else:
                tokens = []
        return loc, tokens

    def _generateDefaultName(self):
        inner = str(self.expr)
        # strip off redundant inner {}'s
        while len(inner) > 1 and inner[0 :: len(inner) - 1] == \"{}\":
            inner = inner[1:-1]
        return \"[\" + inner + \"]\"


Optional = Opt


class SkipTo(ParseElementEnhance):
    \"\"\"
    Token for skipping over all undefined text until the matched
    expression is found.

    Parameters:
    - ``expr`` - target expression marking the end of the data to be skipped
    - ``include`` - if ``True``, the target expression is also parsed
      (the skipped text and target expression are returned as a 2-element
      list) (default= ``False``).
    - ``ignore`` - (default= ``None``) used to define grammars (typically quoted strings and
      comments) that might contain false matches to the target expression
    - ``fail_on`` - (default= ``None``) define expressions that are not allowed to be
      included in the skipped test; if found before the target expression is found,
      the :class:`SkipTo` is not a match

    Example::

        report = '''
            Outstanding Issues Report - 1 Jan 2000

               # | Severity | Description                               |  Days Open
            -----+----------+-------------------------------------------+-----------
             101 | Critical | Intermittent system crash                 |          6
              94 | Cosmetic | Spelling error on Login ('log|n')         |         14
              79 | Minor    | System slow when running too many reports |         47
            '''
        integer = Word(nums)
        SEP = Suppress('|')
        # use SkipTo to simply match everything up until the next SEP
        # - ignore quoted strings, so that a '|' character inside a quoted string does not match
        # - parse action will call token.strip() for each matched token, i.e., the description body
        string_data = SkipTo(SEP, ignore=quoted_string)
        string_data.set_parse_action(token_map(str.strip))
        ticket_expr = (integer(\"issue_num\") + SEP
                      + string_data(\"sev\") + SEP
                      + string_data(\"desc\") + SEP
                      + integer(\"days_open\"))

        for tkt in ticket_expr.search_string(report):
            print tkt.dump()

    prints::

        ['101', 'Critical', 'Intermittent system crash', '6']
        - days_open: '6'
        - desc: 'Intermittent system crash'
        - issue_num: '101'
        - sev: 'Critical'
        ['94', 'Cosmetic', \"Spelling error on Login ('log|n')\", '14']
        - days_open: '14'
        - desc: \"Spelling error on Login ('log|n')\"
        - issue_num: '94'
        - sev: 'Cosmetic'
        ['79', 'Minor', 'System slow when running too many reports', '47']
        - days_open: '47'
        - desc: 'System slow when running too many reports'
        - issue_num: '79'
        - sev: 'Minor'
    \"\"\"

    def __init__(
        self,
        other: Union[ParserElement, str],
        include: bool = False,
        ignore: bool = None,
        fail_on: typing.Optional[Union[ParserElement, str]] = None,
        *,
        failOn: Union[ParserElement, str] = None,
    ):
        super().__init__(other)
        failOn = failOn or fail_on
        self.ignoreExpr = ignore
        self.mayReturnEmpty = True
        self.mayIndexError = False
        self.includeMatch = include
        self.saveAsList = False
        if isinstance(failOn, str_type):
            self.failOn = self._literalStringClass(failOn)
        else:
            self.failOn = failOn
        self.errmsg = \"No match found for \" + str(self.expr)

    def parseImpl(self, instring, loc, doActions=True):
        startloc = loc
        instrlen = len(instring)
        self_expr_parse = self.expr._parse
        self_failOn_canParseNext = (
            self.failOn.canParseNext if self.failOn is not None else None
        )
        self_ignoreExpr_tryParse = (
            self.ignoreExpr.tryParse if self.ignoreExpr is not None else None
        )

        tmploc = loc
        while tmploc <= instrlen:
            if self_failOn_canParseNext is not None:
                # break if failOn expression matches
                if self_failOn_canParseNext(instring, tmploc):
                    break

            if self_ignoreExpr_tryParse is not None:
                # advance past ignore expressions
                while 1:
                    try:
                        tmploc = self_ignoreExpr_tryParse(instring, tmploc)
                    except ParseBaseException:
                        break

            try:
                self_expr_parse(instring, tmploc, doActions=False, callPreParse=False)
            except (ParseException, IndexError):
                # no match, advance loc in string
                tmploc += 1
            else:
                # matched skipto expr, done
                break

        else:
            # ran off the end of the input string without matching skipto expr, fail
            raise ParseException(instring, loc, self.errmsg, self)

        # build up return values
        loc = tmploc
        skiptext = instring[startloc:loc]
        skipresult = ParseResults(skiptext)

        if self.includeMatch:
            loc, mat = self_expr_parse(instring, loc, doActions, callPreParse=False)
            skipresult += mat

        return loc, skipresult


class Forward(ParseElementEnhance):
    \"\"\"
    Forward declaration of an expression to be defined later -
    used for recursive grammars, such as algebraic infix notation.
    When the expression is known, it is assigned to the ``Forward``
    variable using the ``'<<'`` operator.

    Note: take care when assigning to ``Forward`` not to overlook
    precedence of operators.

    Specifically, ``'|'`` has a lower precedence than ``'<<'``, so that::

        fwd_expr << a | b | c

    will actually be evaluated as::

        (fwd_expr << a) | b | c

    thereby leaving b and c out as parseable alternatives.  It is recommended that you
    explicitly group the values inserted into the ``Forward``::

        fwd_expr << (a | b | c)

    Converting to use the ``'<<='`` operator instead will avoid this problem.

    See :class:`ParseResults.pprint` for an example of a recursive
    parser created using ``Forward``.
    \"\"\"

    def __init__(self, other: typing.Optional[Union[ParserElement, str]] = None):
        self.caller_frame = traceback.extract_stack(limit=2)[0]
        super().__init__(other, savelist=False)
        self.lshift_line = None

    def __lshift__(self, other):
        if hasattr(self, \"caller_frame\"):
            del self.caller_frame
        if isinstance(other, str_type):
            other = self._literalStringClass(other)
        self.expr = other
        self.mayIndexError = self.expr.mayIndexError
        self.mayReturnEmpty = self.expr.mayReturnEmpty
        self.set_whitespace_chars(
            self.expr.whiteChars, copy_defaults=self.expr.copyDefaultWhiteChars
        )
        self.skipWhitespace = self.expr.skipWhitespace
        self.saveAsList = self.expr.saveAsList
        self.ignoreExprs.extend(self.expr.ignoreExprs)
        self.lshift_line = traceback.extract_stack(limit=2)[-2]
        return self

    def __ilshift__(self, other):
        return self << other

    def __or__(self, other):
        caller_line = traceback.extract_stack(limit=2)[-2]
        if (
            __diag__.warn_on_match_first_with_lshift_operator
            and caller_line == self.lshift_line
            and Diagnostics.warn_on_match_first_with_lshift_operator
            not in self.suppress_warnings_
        ):
            warnings.warn(
                \"using '<<' operator with '|' is probably an error, use '<<='\",
                stacklevel=2,
            )
        ret = super().__or__(other)
        return ret

    def __del__(self):
        # see if we are getting dropped because of '=' reassignment of var instead of '<<=' or '<<'
        if (
            self.expr is None
            and __diag__.warn_on_assignment_to_Forward
            and Diagnostics.warn_on_assignment_to_Forward not in self.suppress_warnings_
        ):
            warnings.warn_explicit(
                \"Forward defined here but no expression attached later using '<<=' or '<<'\",
                UserWarning,
                filename=self.caller_frame.filename,
                lineno=self.caller_frame.lineno,
            )

    def parseImpl(self, instring, loc, doActions=True):
        if (
            self.expr is None
            and __diag__.warn_on_parse_using_empty_Forward
            and Diagnostics.warn_on_parse_using_empty_Forward
            not in self.suppress_warnings_
        ):
            # walk stack until parse_string, scan_string, search_string, or transform_string is found
            parse_fns = [
                \"parse_string\",
                \"scan_string\",
                \"search_string\",
                \"transform_string\",
            ]
            tb = traceback.extract_stack(limit=200)
            for i, frm in enumerate(reversed(tb), start=1):
                if frm.name in parse_fns:
                    stacklevel = i + 1
                    break
            else:
                stacklevel = 2
            warnings.warn(
                \"Forward expression was never assigned a value, will not parse any input\",
                stacklevel=stacklevel,
            )
        if not ParserElement._left_recursion_enabled:
            return super().parseImpl(instring, loc, doActions)
        # ## Bounded Recursion algorithm ##
        # Recursion only needs to be processed at ``Forward`` elements, since they are
        # the only ones that can actually refer to themselves. The general idea is
        # to handle recursion stepwise: We start at no recursion, then recurse once,
        # recurse twice, ..., until more recursion offers no benefit (we hit the bound).
        #
        # The \"trick\" here is that each ``Forward`` gets evaluated in two contexts
        # - to *match* a specific recursion level, and
        # - to *search* the bounded recursion level
        # and the two run concurrently. The *search* must *match* each recursion level
        # to find the best possible match. This is handled by a memo table, which
        # provides the previous match to the next level match attempt.
        #
        # See also \"Left Recursion in Parsing Expression Grammars\", Medeiros et al.
        #
        # There is a complication since we not only *parse* but also *transform* via
        # actions: We do not want to run the actions too often while expanding. Thus,
        # we expand using `doActions=False` and only run `doActions=True` if the next
        # recursion level is acceptable.
        with ParserElement.recursion_lock:
            memo = ParserElement.recursion_memos
            try:
                # we are parsing at a specific recursion expansion - use it as-is
                prev_loc, prev_result = memo[loc, self, doActions]
                if isinstance(prev_result, Exception):
                    raise prev_result
                return prev_loc, prev_result.copy()
            except KeyError:
                act_key = (loc, self, True)
                peek_key = (loc, self, False)
                # we are searching for the best recursion expansion - keep on improving
                # both `doActions` cases must be tracked separately here!
                prev_loc, prev_peek = memo[peek_key] = (
                    loc - 1,
                    ParseException(
                        instring, loc, \"Forward recursion without base case\", self
                    ),
                )
                if doActions:
                    memo[act_key] = memo[peek_key]
                while True:
                    try:
                        new_loc, new_peek = super().parseImpl(instring, loc, False)
                    except ParseException:
                        # we failed before getting any match – do not hide the error
                        if isinstance(prev_peek, Exception):
                            raise
                        new_loc, new_peek = prev_loc, prev_peek
                    # the match did not get better: we are done
                    if new_loc <= prev_loc:
                        if doActions:
                            # replace the match for doActions=False as well,
                            # in case the action did backtrack
                            prev_loc, prev_result = memo[peek_key] = memo[act_key]
                            del memo[peek_key], memo[act_key]
                            return prev_loc, prev_result.copy()
                        del memo[peek_key]
                        return prev_loc, prev_peek.copy()
                    # the match did get better: see if we can improve further
                    else:
                        if doActions:
                            try:
                                memo[act_key] = super().parseImpl(instring, loc, True)
                            except ParseException as e:
                                memo[peek_key] = memo[act_key] = (new_loc, e)
                                raise
                        prev_loc, prev_peek = memo[peek_key] = new_loc, new_peek

    def leave_whitespace(self, recursive: bool = True) -> ParserElement:
        self.skipWhitespace = False
        return self

    def ignore_whitespace(self, recursive: bool = True) -> ParserElement:
        self.skipWhitespace = True
        return self

    def streamline(self) -> ParserElement:
        if not self.streamlined:
            self.streamlined = True
            if self.expr is not None:
                self.expr.streamline()
        return self

    def validate(self, validateTrace=None) -> None:
        if validateTrace is None:
            validateTrace = []

        if self not in validateTrace:
            tmp = validateTrace[:] + [self]
            if self.expr is not None:
                self.expr.validate(tmp)
        self._checkRecursion([])

    def _generateDefaultName(self):
        # Avoid infinite recursion by setting a temporary _defaultName
        self._defaultName = \": ...\"

        # Use the string representation of main expression.
        retString = \"...\"
        try:
            if self.expr is not None:
                retString = str(self.expr)[:1000]
            else:
                retString = \"None\"
        finally:
            return self.__class__.__name__ + \": \" + retString

    def copy(self) -> ParserElement:
        if self.expr is not None:
            return super().copy()
        else:
            ret = Forward()
            ret <<= self
            return ret

    def _setResultsName(self, name, list_all_matches=False):
        if (
            __diag__.warn_name_set_on_empty_Forward
            and Diagnostics.warn_name_set_on_empty_Forward
            not in self.suppress_warnings_
        ):
            if self.expr is None:
                warnings.warn(
                    \"{}: setting results name {!r} on {} expression \"
                    \"that has no contained expression\".format(
                        \"warn_name_set_on_empty_Forward\", name, type(self).__name__
                    ),
                    stacklevel=3,
                )

        return super()._setResultsName(name, list_all_matches)

    ignoreWhitespace = ignore_whitespace
    leaveWhitespace = leave_whitespace


class TokenConverter(ParseElementEnhance):
    \"\"\"
    Abstract subclass of :class:`ParseExpression`, for converting parsed results.
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str], savelist=False):
        super().__init__(expr)  # , savelist)
        self.saveAsList = False


class Combine(TokenConverter):
    \"\"\"Converter to concatenate all matching tokens to a single string.
    By default, the matching patterns must also be contiguous in the
    input string; this can be disabled by specifying
    ``'adjacent=False'`` in the constructor.

    Example::

        real = Word(nums) + '.' + Word(nums)
        print(real.parse_string('3.1416')) # -> ['3', '.', '1416']
        # will also erroneously match the following
        print(real.parse_string('3. 1416')) # -> ['3', '.', '1416']

        real = Combine(Word(nums) + '.' + Word(nums))
        print(real.parse_string('3.1416')) # -> ['3.1416']
        # no match when there are internal spaces
        print(real.parse_string('3. 1416')) # -> Exception: Expected W:(0123...)
    \"\"\"

    def __init__(
        self,
        expr: ParserElement,
        join_string: str = \"\",
        adjacent: bool = True,
        *,
        joinString: typing.Optional[str] = None,
    ):
        super().__init__(expr)
        joinString = joinString if joinString is not None else join_string
        # suppress whitespace-stripping in contained parse expressions, but re-enable it on the Combine itself
        if adjacent:
            self.leave_whitespace()
        self.adjacent = adjacent
        self.skipWhitespace = True
        self.joinString = joinString
        self.callPreparse = True

    def ignore(self, other) -> ParserElement:
        if self.adjacent:
            ParserElement.ignore(self, other)
        else:
            super().ignore(other)
        return self

    def postParse(self, instring, loc, tokenlist):
        retToks = tokenlist.copy()
        del retToks[:]
        retToks += ParseResults(
            [\"\".join(tokenlist._asStringList(self.joinString))], modal=self.modalResults
        )

        if self.resultsName and retToks.haskeys():
            return [retToks]
        else:
            return retToks


class Group(TokenConverter):
    \"\"\"Converter to return the matched tokens as a list - useful for
    returning tokens of :class:`ZeroOrMore` and :class:`OneOrMore` expressions.

    The optional ``aslist`` argument when set to True will return the
    parsed tokens as a Python list instead of a pyparsing ParseResults.

    Example::

        ident = Word(alphas)
        num = Word(nums)
        term = ident | num
        func = ident + Opt(delimited_list(term))
        print(func.parse_string(\"fn a, b, 100\"))
        # -> ['fn', 'a', 'b', '100']

        func = ident + Group(Opt(delimited_list(term)))
        print(func.parse_string(\"fn a, b, 100\"))
        # -> ['fn', ['a', 'b', '100']]
    \"\"\"

    def __init__(self, expr: ParserElement, aslist: bool = False):
        super().__init__(expr)
        self.saveAsList = True
        self._asPythonList = aslist

    def postParse(self, instring, loc, tokenlist):
        if self._asPythonList:
            return ParseResults.List(
                tokenlist.asList()
                if isinstance(tokenlist, ParseResults)
                else list(tokenlist)
            )
        else:
            return [tokenlist]


class Dict(TokenConverter):
    \"\"\"Converter to return a repetitive expression as a list, but also
    as a dictionary. Each element can also be referenced using the first
    token in the expression as its key. Useful for tabular report
    scraping when the first column can be used as a item key.

    The optional ``asdict`` argument when set to True will return the
    parsed tokens as a Python dict instead of a pyparsing ParseResults.

    Example::

        data_word = Word(alphas)
        label = data_word + FollowedBy(':')

        text = \"shape: SQUARE posn: upper left color: light blue texture: burlap\"
        attr_expr = (label + Suppress(':') + OneOrMore(data_word, stop_on=label).set_parse_action(' '.join))

        # print attributes as plain groups
        print(attr_expr[1, ...].parse_string(text).dump())

        # instead of OneOrMore(expr), parse using Dict(Group(expr)[1, ...]) - Dict will auto-assign names
        result = Dict(Group(attr_expr)[1, ...]).parse_string(text)
        print(result.dump())

        # access named fields as dict entries, or output as dict
        print(result['shape'])
        print(result.as_dict())

    prints::

        ['shape', 'SQUARE', 'posn', 'upper left', 'color', 'light blue', 'texture', 'burlap']
        [['shape', 'SQUARE'], ['posn', 'upper left'], ['color', 'light blue'], ['texture', 'burlap']]
        - color: 'light blue'
        - posn: 'upper left'
        - shape: 'SQUARE'
        - texture: 'burlap'
        SQUARE
        {'color': 'light blue', 'posn': 'upper left', 'texture': 'burlap', 'shape': 'SQUARE'}

    See more examples at :class:`ParseResults` of accessing fields by results name.
    \"\"\"

    def __init__(self, expr: ParserElement, asdict: bool = False):
        super().__init__(expr)
        self.saveAsList = True
        self._asPythonDict = asdict

    def postParse(self, instring, loc, tokenlist):
        for i, tok in enumerate(tokenlist):
            if len(tok) == 0:
                continue

            ikey = tok[0]
            if isinstance(ikey, int):
                ikey = str(ikey).strip()

            if len(tok) == 1:
                tokenlist[ikey] = _ParseResultsWithOffset(\"\", i)

            elif len(tok) == 2 and not isinstance(tok[1], ParseResults):
                tokenlist[ikey] = _ParseResultsWithOffset(tok[1], i)

            else:
                try:
                    dictvalue = tok.copy()  # ParseResults(i)
                except Exception:
                    exc = TypeError(
                        \"could not extract dict values from parsed results\"
                        \" - Dict expression must contain Grouped expressions\"
                    )
                    raise exc from None

                del dictvalue[0]

                if len(dictvalue) != 1 or (
                    isinstance(dictvalue, ParseResults) and dictvalue.haskeys()
                ):
                    tokenlist[ikey] = _ParseResultsWithOffset(dictvalue, i)
                else:
                    tokenlist[ikey] = _ParseResultsWithOffset(dictvalue[0], i)

        if self._asPythonDict:
            return [tokenlist.as_dict()] if self.resultsName else tokenlist.as_dict()
        else:
            return [tokenlist] if self.resultsName else tokenlist


class Suppress(TokenConverter):
    \"\"\"Converter for ignoring the results of a parsed expression.

    Example::

        source = \"a, b, c,d\"
        wd = Word(alphas)
        wd_list1 = wd + (',' + wd)[...]
        print(wd_list1.parse_string(source))

        # often, delimiters that are useful during parsing are just in the
        # way afterward - use Suppress to keep them out of the parsed output
        wd_list2 = wd + (Suppress(',') + wd)[...]
        print(wd_list2.parse_string(source))

        # Skipped text (using '...') can be suppressed as well
        source = \"lead in START relevant text END trailing text\"
        start_marker = Keyword(\"START\")
        end_marker = Keyword(\"END\")
        find_body = Suppress(...) + start_marker + ... + end_marker
        print(find_body.parse_string(source)

    prints::

        ['a', ',', 'b', ',', 'c', ',', 'd']
        ['a', 'b', 'c', 'd']
        ['START', 'relevant text ', 'END']

    (See also :class:`delimited_list`.)
    \"\"\"

    def __init__(self, expr: Union[ParserElement, str], savelist: bool = False):
        if expr is ...:
            expr = _PendingSkip(NoMatch())
        super().__init__(expr)

    def __add__(self, other) -> \"ParserElement\":
        if isinstance(self.expr, _PendingSkip):
            return Suppress(SkipTo(other)) + other
        else:
            return super().__add__(other)

    def __sub__(self, other) -> \"ParserElement\":
        if isinstance(self.expr, _PendingSkip):
            return Suppress(SkipTo(other)) - other
        else:
            return super().__sub__(other)

    def postParse(self, instring, loc, tokenlist):
        return []

    def suppress(self) -> ParserElement:
        return self


def trace_parse_action(f: ParseAction) -> ParseAction:
    \"\"\"Decorator for debugging parse actions.

    When the parse action is called, this decorator will print
    ``\">> entering method-name(line:<current_source_line>, <parse_location>, <matched_tokens>)\"``.
    When the parse action completes, the decorator will print
    ``\"<<\"`` followed by the returned value, or any exception that the parse action raised.

    Example::

        wd = Word(alphas)

        @trace_parse_action
        def remove_duplicate_chars(tokens):
            return ''.join(sorted(set(''.join(tokens))))

        wds = wd[1, ...].set_parse_action(remove_duplicate_chars)
        print(wds.parse_string(\"slkdjs sld sldd sdlf sdljf\"))

    prints::

        >>entering remove_duplicate_chars(line: 'slkdjs sld sldd sdlf sdljf', 0, (['slkdjs', 'sld', 'sldd', 'sdlf', 'sdljf'], {}))
        <<leaving remove_duplicate_chars (ret: 'dfjkls')
        ['dfjkls']
    \"\"\"
    f = _trim_arity(f)

    def z(*paArgs):
        thisFunc = f.__name__
        s, l, t = paArgs[-3:]
        if len(paArgs) > 3:
            thisFunc = paArgs[0].__class__.__name__ + \".\" + thisFunc
        sys.stderr.write(
            \">>entering {}(line: {!r}, {}, {!r})\\n\".format(thisFunc, line(l, s), l, t)
        )
        try:
            ret = f(*paArgs)
        except Exception as exc:
            sys.stderr.write(\"<<leaving {} (exception: {})\\n\".format(thisFunc, exc))
            raise
        sys.stderr.write(\"<<leaving {} (ret: {!r})\\n\".format(thisFunc, ret))
        return ret

    z.__name__ = f.__name__
    return z


# convenience constants for positional expressions
empty = Empty().set_name(\"empty\")
line_start = LineStart().set_name(\"line_start\")
line_end = LineEnd().set_name(\"line_end\")
string_start = StringStart().set_name(\"string_start\")
string_end = StringEnd().set_name(\"string_end\")

_escapedPunc = Word(_bslash, r\"\\[]-*.$+^?()~ \", exact=2).set_parse_action(
    lambda s, l, t: t[0][1]
)
_escapedHexChar = Regex(r\"\\\\0?[xX][0-9a-fA-F]+\").set_parse_action(
    lambda s, l, t: chr(int(t[0].lstrip(r\"\\0x\"), 16))
)
_escapedOctChar = Regex(r\"\\\\0[0-7]+\").set_parse_action(
    lambda s, l, t: chr(int(t[0][1:], 8))
)
_singleChar = (
    _escapedPunc | _escapedHexChar | _escapedOctChar | CharsNotIn(r\"\\]\", exact=1)
)
_charRange = Group(_singleChar + Suppress(\"-\") + _singleChar)
_reBracketExpr = (
    Literal(\"[\")
    + Opt(\"^\").set_results_name(\"negate\")
    + Group(OneOrMore(_charRange | _singleChar)).set_results_name(\"body\")
    + \"]\"
)


def srange(s: str) -> str:
    r\"\"\"Helper to easily define string ranges for use in :class:`Word`
    construction. Borrows syntax from regexp ``'[]'`` string range
    definitions::

        srange(\"[0-9]\")   -> \"0123456789\"
        srange(\"[a-z]\")   -> \"abcdefghijklmnopqrstuvwxyz\"
        srange(\"[a-z$_]\") -> \"abcdefghijklmnopqrstuvwxyz$_\"

    The input string must be enclosed in []'s, and the returned string
    is the expanded character set joined into a single string. The
    values enclosed in the []'s may be:

    - a single character
    - an escaped character with a leading backslash (such as ``\\-``
      or ``\\]``)
    - an escaped hex character with a leading ``'\\x'``
      (``\\x21``, which is a ``'!'`` character) (``\\0x##``
      is also supported for backwards compatibility)
    - an escaped octal character with a leading ``'\\0'``
      (``\\041``, which is a ``'!'`` character)
    - a range of any of the above, separated by a dash (``'a-z'``,
      etc.)
    - any combination of the above (``'aeiouy'``,
      ``'a-zA-Z0-9_$'``, etc.)
    \"\"\"
    _expanded = (
        lambda p: p
        if not isinstance(p, ParseResults)
        else \"\".join(chr(c) for c in range(ord(p[0]), ord(p[1]) + 1))
    )
    try:
        return \"\".join(_expanded(part) for part in _reBracketExpr.parse_string(s).body)
    except Exception:
        return \"\"


def token_map(func, *args) -> ParseAction:
    \"\"\"Helper to define a parse action by mapping a function to all
    elements of a :class:`ParseResults` list. If any additional args are passed,
    they are forwarded to the given function as additional arguments
    after the token, as in
    ``hex_integer = Word(hexnums).set_parse_action(token_map(int, 16))``,
    which will convert the parsed data to an integer using base 16.

    Example (compare the last to example in :class:`ParserElement.transform_string`::

        hex_ints = Word(hexnums)[1, ...].set_parse_action(token_map(int, 16))
        hex_ints.run_tests('''
            00 11 22 aa FF 0a 0d 1a
            ''')

        upperword = Word(alphas).set_parse_action(token_map(str.upper))
        upperword[1, ...].run_tests('''
            my kingdom for a horse
            ''')

        wd = Word(alphas).set_parse_action(token_map(str.title))
        wd[1, ...].set_parse_action(' '.join).run_tests('''
            now is the winter of our discontent made glorious summer by this sun of york
            ''')

    prints::

        00 11 22 aa FF 0a 0d 1a
        [0, 17, 34, 170, 255, 10, 13, 26]

        my kingdom for a horse
        ['MY', 'KINGDOM', 'FOR', 'A', 'HORSE']

        now is the winter of our discontent made glorious summer by this sun of york
        ['Now Is The Winter Of Our Discontent Made Glorious Summer By This Sun Of York']
    \"\"\"

    def pa(s, l, t):
        return [func(tokn, *args) for tokn in t]

    func_name = getattr(func, \"__name__\", getattr(func, \"__class__\").__name__)
    pa.__name__ = func_name

    return pa


def autoname_elements() -> None:
    \"\"\"
    Utility to simplify mass-naming of parser elements, for
    generating railroad diagram with named subdiagrams.
    \"\"\"
    for name, var in sys._getframe().f_back.f_locals.items():
        if isinstance(var, ParserElement) and not var.customName:
            var.set_name(name)


dbl_quoted_string = Combine(
    Regex(r'\"(?:[^\"\\n\\r\\\\]|(?:\"\")|(?:\\\\(?:[^x]|x[0-9a-fA-F]+)))*') + '\"'
).set_name(\"string enclosed in double quotes\")

sgl_quoted_string = Combine(
    Regex(r\"'(?:[^'\\n\\r\\\\]|(?:'')|(?:\\\\(?:[^x]|x[0-9a-fA-F]+)))*\") + \"'\"
).set_name(\"string enclosed in single quotes\")

quoted_string = Combine(
    Regex(r'\"(?:[^\"\\n\\r\\\\]|(?:\"\")|(?:\\\\(?:[^x]|x[0-9a-fA-F]+)))*') + '\"'
    | Regex(r\"'(?:[^'\\n\\r\\\\]|(?:'')|(?:\\\\(?:[^x]|x[0-9a-fA-F]+)))*\") + \"'\"
).set_name(\"quotedString using single or double quotes\")

unicode_string = Combine(\"u\" + quoted_string.copy()).set_name(\"unicode string literal\")


alphas8bit = srange(r\"[\\0xc0-\\0xd6\\0xd8-\\0xf6\\0xf8-\\0xff]\")
punc8bit = srange(r\"[\\0xa1-\\0xbf\\0xd7\\0xf7]\")

# build list of built-in expressions, for future reference if a global default value
# gets updated
_builtin_exprs: List[ParserElement] = [
    v for v in vars().values() if isinstance(v, ParserElement)
]

# backward compatibility names
tokenMap = token_map
conditionAsParseAction = condition_as_parse_action
nullDebugAction = null_debug_action
sglQuotedString = sgl_quoted_string
dblQuotedString = dbl_quoted_string
quotedString = quoted_string
unicodeString = unicode_string
lineStart = line_start
lineEnd = line_end
stringStart = string_start
stringEnd = string_end
traceParseAction = trace_parse_action

"""
module_dict["pyparsing/helpers.py"]="""
# helpers.py
import html.entities
import re
import typing

from . import __diag__
from .core import *
from .util import _bslash, _flatten, _escape_regex_range_chars


#
# global helpers
#
def delimited_list(
    expr: Union[str, ParserElement],
    delim: Union[str, ParserElement] = \",\",
    combine: bool = False,
    min: typing.Optional[int] = None,
    max: typing.Optional[int] = None,
    *,
    allow_trailing_delim: bool = False,
) -> ParserElement:
    \"\"\"Helper to define a delimited list of expressions - the delimiter
    defaults to ','. By default, the list elements and delimiters can
    have intervening whitespace, and comments, but this can be
    overridden by passing ``combine=True`` in the constructor. If
    ``combine`` is set to ``True``, the matching tokens are
    returned as a single token string, with the delimiters included;
    otherwise, the matching tokens are returned as a list of tokens,
    with the delimiters suppressed.

    If ``allow_trailing_delim`` is set to True, then the list may end with
    a delimiter.

    Example::

        delimited_list(Word(alphas)).parse_string(\"aa,bb,cc\") # -> ['aa', 'bb', 'cc']
        delimited_list(Word(hexnums), delim=':', combine=True).parse_string(\"AA:BB:CC:DD:EE\") # -> ['AA:BB:CC:DD:EE']
    \"\"\"
    if isinstance(expr, str_type):
        expr = ParserElement._literalStringClass(expr)

    dlName = \"{expr} [{delim} {expr}]...{end}\".format(
        expr=str(expr.copy().streamline()),
        delim=str(delim),
        end=\" [{}]\".format(str(delim)) if allow_trailing_delim else \"\",
    )

    if not combine:
        delim = Suppress(delim)

    if min is not None:
        if min < 1:
            raise ValueError(\"min must be greater than 0\")
        min -= 1
    if max is not None:
        if min is not None and max <= min:
            raise ValueError(\"max must be greater than, or equal to min\")
        max -= 1
    delimited_list_expr = expr + (delim + expr)[min, max]

    if allow_trailing_delim:
        delimited_list_expr += Opt(delim)

    if combine:
        return Combine(delimited_list_expr).set_name(dlName)
    else:
        return delimited_list_expr.set_name(dlName)


def counted_array(
    expr: ParserElement,
    int_expr: typing.Optional[ParserElement] = None,
    *,
    intExpr: typing.Optional[ParserElement] = None,
) -> ParserElement:
    \"\"\"Helper to define a counted list of expressions.

    This helper defines a pattern of the form::

        integer expr expr expr...

    where the leading integer tells how many expr expressions follow.
    The matched tokens returns the array of expr tokens as a list - the
    leading count token is suppressed.

    If ``int_expr`` is specified, it should be a pyparsing expression
    that produces an integer value.

    Example::

        counted_array(Word(alphas)).parse_string('2 ab cd ef')  # -> ['ab', 'cd']

        # in this parser, the leading integer value is given in binary,
        # '10' indicating that 2 values are in the array
        binary_constant = Word('01').set_parse_action(lambda t: int(t[0], 2))
        counted_array(Word(alphas), int_expr=binary_constant).parse_string('10 ab cd ef')  # -> ['ab', 'cd']

        # if other fields must be parsed after the count but before the
        # list items, give the fields results names and they will
        # be preserved in the returned ParseResults:
        count_with_metadata = integer + Word(alphas)(\"type\")
        typed_array = counted_array(Word(alphanums), int_expr=count_with_metadata)(\"items\")
        result = typed_array.parse_string(\"3 bool True True False\")
        print(result.dump())

        # prints
        # ['True', 'True', 'False']
        # - items: ['True', 'True', 'False']
        # - type: 'bool'
    \"\"\"
    intExpr = intExpr or int_expr
    array_expr = Forward()

    def count_field_parse_action(s, l, t):
        nonlocal array_expr
        n = t[0]
        array_expr <<= (expr * n) if n else Empty()
        # clear list contents, but keep any named results
        del t[:]

    if intExpr is None:
        intExpr = Word(nums).set_parse_action(lambda t: int(t[0]))
    else:
        intExpr = intExpr.copy()
    intExpr.set_name(\"arrayLen\")
    intExpr.add_parse_action(count_field_parse_action, call_during_try=True)
    return (intExpr + array_expr).set_name(\"(len) \" + str(expr) + \"...\")


def match_previous_literal(expr: ParserElement) -> ParserElement:
    \"\"\"Helper to define an expression that is indirectly defined from
    the tokens matched in a previous expression, that is, it looks for
    a 'repeat' of a previous expression.  For example::

        first = Word(nums)
        second = match_previous_literal(first)
        match_expr = first + \":\" + second

    will match ``\"1:1\"``, but not ``\"1:2\"``.  Because this
    matches a previous literal, will also match the leading
    ``\"1:1\"`` in ``\"1:10\"``. If this is not desired, use
    :class:`match_previous_expr`. Do *not* use with packrat parsing
    enabled.
    \"\"\"
    rep = Forward()

    def copy_token_to_repeater(s, l, t):
        if t:
            if len(t) == 1:
                rep << t[0]
            else:
                # flatten t tokens
                tflat = _flatten(t.as_list())
                rep << And(Literal(tt) for tt in tflat)
        else:
            rep << Empty()

    expr.add_parse_action(copy_token_to_repeater, callDuringTry=True)
    rep.set_name(\"(prev) \" + str(expr))
    return rep


def match_previous_expr(expr: ParserElement) -> ParserElement:
    \"\"\"Helper to define an expression that is indirectly defined from
    the tokens matched in a previous expression, that is, it looks for
    a 'repeat' of a previous expression.  For example::

        first = Word(nums)
        second = match_previous_expr(first)
        match_expr = first + \":\" + second

    will match ``\"1:1\"``, but not ``\"1:2\"``.  Because this
    matches by expressions, will *not* match the leading ``\"1:1\"``
    in ``\"1:10\"``; the expressions are evaluated first, and then
    compared, so ``\"1\"`` is compared with ``\"10\"``. Do *not* use
    with packrat parsing enabled.
    \"\"\"
    rep = Forward()
    e2 = expr.copy()
    rep <<= e2

    def copy_token_to_repeater(s, l, t):
        matchTokens = _flatten(t.as_list())

        def must_match_these_tokens(s, l, t):
            theseTokens = _flatten(t.as_list())
            if theseTokens != matchTokens:
                raise ParseException(
                    s, l, \"Expected {}, found{}\".format(matchTokens, theseTokens)
                )

        rep.set_parse_action(must_match_these_tokens, callDuringTry=True)

    expr.add_parse_action(copy_token_to_repeater, callDuringTry=True)
    rep.set_name(\"(prev) \" + str(expr))
    return rep


def one_of(
    strs: Union[typing.Iterable[str], str],
    caseless: bool = False,
    use_regex: bool = True,
    as_keyword: bool = False,
    *,
    useRegex: bool = True,
    asKeyword: bool = False,
) -> ParserElement:
    \"\"\"Helper to quickly define a set of alternative :class:`Literal` s,
    and makes sure to do longest-first testing when there is a conflict,
    regardless of the input order, but returns
    a :class:`MatchFirst` for best performance.

    Parameters:

    - ``strs`` - a string of space-delimited literals, or a collection of
      string literals
    - ``caseless`` - treat all literals as caseless - (default= ``False``)
    - ``use_regex`` - as an optimization, will
      generate a :class:`Regex` object; otherwise, will generate
      a :class:`MatchFirst` object (if ``caseless=True`` or ``asKeyword=True``, or if
      creating a :class:`Regex` raises an exception) - (default= ``True``)
    - ``as_keyword`` - enforce :class:`Keyword`-style matching on the
      generated expressions - (default= ``False``)
    - ``asKeyword`` and ``useRegex`` are retained for pre-PEP8 compatibility,
      but will be removed in a future release

    Example::

        comp_oper = one_of(\"< = > <= >= !=\")
        var = Word(alphas)
        number = Word(nums)
        term = var | number
        comparison_expr = term + comp_oper + term
        print(comparison_expr.search_string(\"B = 12  AA=23 B<=AA AA>12\"))

    prints::

        [['B', '=', '12'], ['AA', '=', '23'], ['B', '<=', 'AA'], ['AA', '>', '12']]
    \"\"\"
    asKeyword = asKeyword or as_keyword
    useRegex = useRegex and use_regex

    if (
        isinstance(caseless, str_type)
        and __diag__.warn_on_multiple_string_args_to_oneof
    ):
        warnings.warn(
            \"More than one string argument passed to one_of, pass\"
            \" choices as a list or space-delimited string\",
            stacklevel=2,
        )

    if caseless:
        isequal = lambda a, b: a.upper() == b.upper()
        masks = lambda a, b: b.upper().startswith(a.upper())
        parseElementClass = CaselessKeyword if asKeyword else CaselessLiteral
    else:
        isequal = lambda a, b: a == b
        masks = lambda a, b: b.startswith(a)
        parseElementClass = Keyword if asKeyword else Literal

    symbols: List[str] = []
    if isinstance(strs, str_type):
        symbols = strs.split()
    elif isinstance(strs, Iterable):
        symbols = list(strs)
    else:
        raise TypeError(\"Invalid argument to one_of, expected string or iterable\")
    if not symbols:
        return NoMatch()

    # reorder given symbols to take care to avoid masking longer choices with shorter ones
    # (but only if the given symbols are not just single characters)
    if any(len(sym) > 1 for sym in symbols):
        i = 0
        while i < len(symbols) - 1:
            cur = symbols[i]
            for j, other in enumerate(symbols[i + 1 :]):
                if isequal(other, cur):
                    del symbols[i + j + 1]
                    break
                elif masks(cur, other):
                    del symbols[i + j + 1]
                    symbols.insert(i, other)
                    break
            else:
                i += 1

    if useRegex:
        re_flags: int = re.IGNORECASE if caseless else 0

        try:
            if all(len(sym) == 1 for sym in symbols):
                # symbols are just single characters, create range regex pattern
                patt = \"[{}]\".format(
                    \"\".join(_escape_regex_range_chars(sym) for sym in symbols)
                )
            else:
                patt = \"|\".join(re.escape(sym) for sym in symbols)

            # wrap with \\b word break markers if defining as keywords
            if asKeyword:
                patt = r\"\\b(?:{})\\b\".format(patt)

            ret = Regex(patt, flags=re_flags).set_name(\" | \".join(symbols))

            if caseless:
                # add parse action to return symbols as specified, not in random
                # casing as found in input string
                symbol_map = {sym.lower(): sym for sym in symbols}
                ret.add_parse_action(lambda s, l, t: symbol_map[t[0].lower()])

            return ret

        except re.error:
            warnings.warn(
                \"Exception creating Regex for one_of, building MatchFirst\", stacklevel=2
            )

    # last resort, just use MatchFirst
    return MatchFirst(parseElementClass(sym) for sym in symbols).set_name(
        \" | \".join(symbols)
    )


def dict_of(key: ParserElement, value: ParserElement) -> ParserElement:
    \"\"\"Helper to easily and clearly define a dictionary by specifying
    the respective patterns for the key and value.  Takes care of
    defining the :class:`Dict`, :class:`ZeroOrMore`, and
    :class:`Group` tokens in the proper order.  The key pattern
    can include delimiting markers or punctuation, as long as they are
    suppressed, thereby leaving the significant key text.  The value
    pattern can include named results, so that the :class:`Dict` results
    can include named token fields.

    Example::

        text = \"shape: SQUARE posn: upper left color: light blue texture: burlap\"
        attr_expr = (label + Suppress(':') + OneOrMore(data_word, stop_on=label).set_parse_action(' '.join))
        print(attr_expr[1, ...].parse_string(text).dump())

        attr_label = label
        attr_value = Suppress(':') + OneOrMore(data_word, stop_on=label).set_parse_action(' '.join)

        # similar to Dict, but simpler call format
        result = dict_of(attr_label, attr_value).parse_string(text)
        print(result.dump())
        print(result['shape'])
        print(result.shape)  # object attribute access works too
        print(result.as_dict())

    prints::

        [['shape', 'SQUARE'], ['posn', 'upper left'], ['color', 'light blue'], ['texture', 'burlap']]
        - color: 'light blue'
        - posn: 'upper left'
        - shape: 'SQUARE'
        - texture: 'burlap'
        SQUARE
        SQUARE
        {'color': 'light blue', 'shape': 'SQUARE', 'posn': 'upper left', 'texture': 'burlap'}
    \"\"\"
    return Dict(OneOrMore(Group(key + value)))


def original_text_for(
    expr: ParserElement, as_string: bool = True, *, asString: bool = True
) -> ParserElement:
    \"\"\"Helper to return the original, untokenized text for a given
    expression.  Useful to restore the parsed fields of an HTML start
    tag into the raw tag text itself, or to revert separate tokens with
    intervening whitespace back to the original matching input text. By
    default, returns astring containing the original parsed text.

    If the optional ``as_string`` argument is passed as
    ``False``, then the return value is
    a :class:`ParseResults` containing any results names that
    were originally matched, and a single token containing the original
    matched text from the input string.  So if the expression passed to
    :class:`original_text_for` contains expressions with defined
    results names, you must set ``as_string`` to ``False`` if you
    want to preserve those results name values.

    The ``asString`` pre-PEP8 argument is retained for compatibility,
    but will be removed in a future release.

    Example::

        src = \"this is test <b> bold <i>text</i> </b> normal text \"
        for tag in (\"b\", \"i\"):
            opener, closer = make_html_tags(tag)
            patt = original_text_for(opener + SkipTo(closer) + closer)
            print(patt.search_string(src)[0])

    prints::

        ['<b> bold <i>text</i> </b>']
        ['<i>text</i>']
    \"\"\"
    asString = asString and as_string

    locMarker = Empty().set_parse_action(lambda s, loc, t: loc)
    endlocMarker = locMarker.copy()
    endlocMarker.callPreparse = False
    matchExpr = locMarker(\"_original_start\") + expr + endlocMarker(\"_original_end\")
    if asString:
        extractText = lambda s, l, t: s[t._original_start : t._original_end]
    else:

        def extractText(s, l, t):
            t[:] = [s[t.pop(\"_original_start\") : t.pop(\"_original_end\")]]

    matchExpr.set_parse_action(extractText)
    matchExpr.ignoreExprs = expr.ignoreExprs
    matchExpr.suppress_warning(Diagnostics.warn_ungrouped_named_tokens_in_collection)
    return matchExpr


def ungroup(expr: ParserElement) -> ParserElement:
    \"\"\"Helper to undo pyparsing's default grouping of And expressions,
    even if all but one are non-empty.
    \"\"\"
    return TokenConverter(expr).add_parse_action(lambda t: t[0])


def locatedExpr(expr: ParserElement) -> ParserElement:
    \"\"\"
    (DEPRECATED - future code should use the Located class)
    Helper to decorate a returned token with its starting and ending
    locations in the input string.

    This helper adds the following results names:

    - ``locn_start`` - location where matched expression begins
    - ``locn_end`` - location where matched expression ends
    - ``value`` - the actual parsed results

    Be careful if the input text contains ``<TAB>`` characters, you
    may want to call :class:`ParserElement.parseWithTabs`

    Example::

        wd = Word(alphas)
        for match in locatedExpr(wd).searchString(\"ljsdf123lksdjjf123lkkjj1222\"):
            print(match)

    prints::

        [[0, 'ljsdf', 5]]
        [[8, 'lksdjjf', 15]]
        [[18, 'lkkjj', 23]]
    \"\"\"
    locator = Empty().set_parse_action(lambda ss, ll, tt: ll)
    return Group(
        locator(\"locn_start\")
        + expr(\"value\")
        + locator.copy().leaveWhitespace()(\"locn_end\")
    )


def nested_expr(
    opener: Union[str, ParserElement] = \"(\",
    closer: Union[str, ParserElement] = \")\",
    content: typing.Optional[ParserElement] = None,
    ignore_expr: ParserElement = quoted_string(),
    *,
    ignoreExpr: ParserElement = quoted_string(),
) -> ParserElement:
    \"\"\"Helper method for defining nested lists enclosed in opening and
    closing delimiters (``\"(\"`` and ``\")\"`` are the default).

    Parameters:
    - ``opener`` - opening character for a nested list
      (default= ``\"(\"``); can also be a pyparsing expression
    - ``closer`` - closing character for a nested list
      (default= ``\")\"``); can also be a pyparsing expression
    - ``content`` - expression for items within the nested lists
      (default= ``None``)
    - ``ignore_expr`` - expression for ignoring opening and closing delimiters
      (default= :class:`quoted_string`)
    - ``ignoreExpr`` - this pre-PEP8 argument is retained for compatibility
      but will be removed in a future release

    If an expression is not provided for the content argument, the
    nested expression will capture all whitespace-delimited content
    between delimiters as a list of separate values.

    Use the ``ignore_expr`` argument to define expressions that may
    contain opening or closing characters that should not be treated as
    opening or closing characters for nesting, such as quoted_string or
    a comment expression.  Specify multiple expressions using an
    :class:`Or` or :class:`MatchFirst`. The default is
    :class:`quoted_string`, but if no expressions are to be ignored, then
    pass ``None`` for this argument.

    Example::

        data_type = one_of(\"void int short long char float double\")
        decl_data_type = Combine(data_type + Opt(Word('*')))
        ident = Word(alphas+'_', alphanums+'_')
        number = pyparsing_common.number
        arg = Group(decl_data_type + ident)
        LPAR, RPAR = map(Suppress, \"()\")

        code_body = nested_expr('{', '}', ignore_expr=(quoted_string | c_style_comment))

        c_function = (decl_data_type(\"type\")
                      + ident(\"name\")
                      + LPAR + Opt(delimited_list(arg), [])(\"args\") + RPAR
                      + code_body(\"body\"))
        c_function.ignore(c_style_comment)

        source_code = '''
            int is_odd(int x) {
                return (x%2);
            }

            int dec_to_hex(char hchar) {
                if (hchar >= '0' && hchar <= '9') {
                    return (ord(hchar)-ord('0'));
                } else {
                    return (10+ord(hchar)-ord('A'));
                }
            }
        '''
        for func in c_function.search_string(source_code):
            print(\"%(name)s (%(type)s) args: %(args)s\" % func)


    prints::

        is_odd (int) args: [['int', 'x']]
        dec_to_hex (int) args: [['char', 'hchar']]
    \"\"\"
    if ignoreExpr != ignore_expr:
        ignoreExpr = ignore_expr if ignoreExpr == quoted_string() else ignoreExpr
    if opener == closer:
        raise ValueError(\"opening and closing strings cannot be the same\")
    if content is None:
        if isinstance(opener, str_type) and isinstance(closer, str_type):
            if len(opener) == 1 and len(closer) == 1:
                if ignoreExpr is not None:
                    content = Combine(
                        OneOrMore(
                            ~ignoreExpr
                            + CharsNotIn(
                                opener + closer + ParserElement.DEFAULT_WHITE_CHARS,
                                exact=1,
                            )
                        )
                    ).set_parse_action(lambda t: t[0].strip())
                else:
                    content = empty.copy() + CharsNotIn(
                        opener + closer + ParserElement.DEFAULT_WHITE_CHARS
                    ).set_parse_action(lambda t: t[0].strip())
            else:
                if ignoreExpr is not None:
                    content = Combine(
                        OneOrMore(
                            ~ignoreExpr
                            + ~Literal(opener)
                            + ~Literal(closer)
                            + CharsNotIn(ParserElement.DEFAULT_WHITE_CHARS, exact=1)
                        )
                    ).set_parse_action(lambda t: t[0].strip())
                else:
                    content = Combine(
                        OneOrMore(
                            ~Literal(opener)
                            + ~Literal(closer)
                            + CharsNotIn(ParserElement.DEFAULT_WHITE_CHARS, exact=1)
                        )
                    ).set_parse_action(lambda t: t[0].strip())
        else:
            raise ValueError(
                \"opening and closing arguments must be strings if no content expression is given\"
            )
    ret = Forward()
    if ignoreExpr is not None:
        ret <<= Group(
            Suppress(opener) + ZeroOrMore(ignoreExpr | ret | content) + Suppress(closer)
        )
    else:
        ret <<= Group(Suppress(opener) + ZeroOrMore(ret | content) + Suppress(closer))
    ret.set_name(\"nested %s%s expression\" % (opener, closer))
    return ret


def _makeTags(tagStr, xml, suppress_LT=Suppress(\"<\"), suppress_GT=Suppress(\">\")):
    \"\"\"Internal helper to construct opening and closing tag expressions, given a tag name\"\"\"
    if isinstance(tagStr, str_type):
        resname = tagStr
        tagStr = Keyword(tagStr, caseless=not xml)
    else:
        resname = tagStr.name

    tagAttrName = Word(alphas, alphanums + \"_-:\")
    if xml:
        tagAttrValue = dbl_quoted_string.copy().set_parse_action(remove_quotes)
        openTag = (
            suppress_LT
            + tagStr(\"tag\")
            + Dict(ZeroOrMore(Group(tagAttrName + Suppress(\"=\") + tagAttrValue)))
            + Opt(\"/\", default=[False])(\"empty\").set_parse_action(
                lambda s, l, t: t[0] == \"/\"
            )
            + suppress_GT
        )
    else:
        tagAttrValue = quoted_string.copy().set_parse_action(remove_quotes) | Word(
            printables, exclude_chars=\">\"
        )
        openTag = (
            suppress_LT
            + tagStr(\"tag\")
            + Dict(
                ZeroOrMore(
                    Group(
                        tagAttrName.set_parse_action(lambda t: t[0].lower())
                        + Opt(Suppress(\"=\") + tagAttrValue)
                    )
                )
            )
            + Opt(\"/\", default=[False])(\"empty\").set_parse_action(
                lambda s, l, t: t[0] == \"/\"
            )
            + suppress_GT
        )
    closeTag = Combine(Literal(\"</\") + tagStr + \">\", adjacent=False)

    openTag.set_name(\"<%s>\" % resname)
    # add start<tagname> results name in parse action now that ungrouped names are not reported at two levels
    openTag.add_parse_action(
        lambda t: t.__setitem__(
            \"start\" + \"\".join(resname.replace(\":\", \" \").title().split()), t.copy()
        )
    )
    closeTag = closeTag(
        \"end\" + \"\".join(resname.replace(\":\", \" \").title().split())
    ).set_name(\"</%s>\" % resname)
    openTag.tag = resname
    closeTag.tag = resname
    openTag.tag_body = SkipTo(closeTag())
    return openTag, closeTag


def make_html_tags(
    tag_str: Union[str, ParserElement]
) -> Tuple[ParserElement, ParserElement]:
    \"\"\"Helper to construct opening and closing tag expressions for HTML,
    given a tag name. Matches tags in either upper or lower case,
    attributes with namespaces and with quoted or unquoted values.

    Example::

        text = '<td>More info at the <a href=\"https://github.com/pyparsing/pyparsing/wiki\">pyparsing</a> wiki page</td>'
        # make_html_tags returns pyparsing expressions for the opening and
        # closing tags as a 2-tuple
        a, a_end = make_html_tags(\"A\")
        link_expr = a + SkipTo(a_end)(\"link_text\") + a_end

        for link in link_expr.search_string(text):
            # attributes in the <A> tag (like \"href\" shown here) are
            # also accessible as named results
            print(link.link_text, '->', link.href)

    prints::

        pyparsing -> https://github.com/pyparsing/pyparsing/wiki
    \"\"\"
    return _makeTags(tag_str, False)


def make_xml_tags(
    tag_str: Union[str, ParserElement]
) -> Tuple[ParserElement, ParserElement]:
    \"\"\"Helper to construct opening and closing tag expressions for XML,
    given a tag name. Matches tags only in the given upper/lower case.

    Example: similar to :class:`make_html_tags`
    \"\"\"
    return _makeTags(tag_str, True)


any_open_tag: ParserElement
any_close_tag: ParserElement
any_open_tag, any_close_tag = make_html_tags(
    Word(alphas, alphanums + \"_:\").set_name(\"any tag\")
)

_htmlEntityMap = {k.rstrip(\";\"): v for k, v in html.entities.html5.items()}
common_html_entity = Regex(\"&(?P<entity>\" + \"|\".join(_htmlEntityMap) + \");\").set_name(
    \"common HTML entity\"
)


def replace_html_entity(t):
    \"\"\"Helper parser action to replace common HTML entities with their special characters\"\"\"
    return _htmlEntityMap.get(t.entity)


class OpAssoc(Enum):
    LEFT = 1
    RIGHT = 2


InfixNotationOperatorArgType = Union[
    ParserElement, str, Tuple[Union[ParserElement, str], Union[ParserElement, str]]
]
InfixNotationOperatorSpec = Union[
    Tuple[
        InfixNotationOperatorArgType,
        int,
        OpAssoc,
        typing.Optional[ParseAction],
    ],
    Tuple[
        InfixNotationOperatorArgType,
        int,
        OpAssoc,
    ],
]


def infix_notation(
    base_expr: ParserElement,
    op_list: List[InfixNotationOperatorSpec],
    lpar: Union[str, ParserElement] = Suppress(\"(\"),
    rpar: Union[str, ParserElement] = Suppress(\")\"),
) -> ParserElement:
    \"\"\"Helper method for constructing grammars of expressions made up of
    operators working in a precedence hierarchy.  Operators may be unary
    or binary, left- or right-associative.  Parse actions can also be
    attached to operator expressions. The generated parser will also
    recognize the use of parentheses to override operator precedences
    (see example below).

    Note: if you define a deep operator list, you may see performance
    issues when using infix_notation. See
    :class:`ParserElement.enable_packrat` for a mechanism to potentially
    improve your parser performance.

    Parameters:
    - ``base_expr`` - expression representing the most basic operand to
      be used in the expression
    - ``op_list`` - list of tuples, one for each operator precedence level
      in the expression grammar; each tuple is of the form ``(op_expr,
      num_operands, right_left_assoc, (optional)parse_action)``, where:

      - ``op_expr`` is the pyparsing expression for the operator; may also
        be a string, which will be converted to a Literal; if ``num_operands``
        is 3, ``op_expr`` is a tuple of two expressions, for the two
        operators separating the 3 terms
      - ``num_operands`` is the number of terms for this operator (must be 1,
        2, or 3)
      - ``right_left_assoc`` is the indicator whether the operator is right
        or left associative, using the pyparsing-defined constants
        ``OpAssoc.RIGHT`` and ``OpAssoc.LEFT``.
      - ``parse_action`` is the parse action to be associated with
        expressions matching this operator expression (the parse action
        tuple member may be omitted); if the parse action is passed
        a tuple or list of functions, this is equivalent to calling
        ``set_parse_action(*fn)``
        (:class:`ParserElement.set_parse_action`)
    - ``lpar`` - expression for matching left-parentheses; if passed as a
      str, then will be parsed as Suppress(lpar). If lpar is passed as
      an expression (such as ``Literal('(')``), then it will be kept in
      the parsed results, and grouped with them. (default= ``Suppress('(')``)
    - ``rpar`` - expression for matching right-parentheses; if passed as a
      str, then will be parsed as Suppress(rpar). If rpar is passed as
      an expression (such as ``Literal(')')``), then it will be kept in
      the parsed results, and grouped with them. (default= ``Suppress(')')``)

    Example::

        # simple example of four-function arithmetic with ints and
        # variable names
        integer = pyparsing_common.signed_integer
        varname = pyparsing_common.identifier

        arith_expr = infix_notation(integer | varname,
            [
            ('-', 1, OpAssoc.RIGHT),
            (one_of('* /'), 2, OpAssoc.LEFT),
            (one_of('+ -'), 2, OpAssoc.LEFT),
            ])

        arith_expr.run_tests('''
            5+3*6
            (5+3)*6
            -2--11
            ''', full_dump=False)

    prints::

        5+3*6
        [[5, '+', [3, '*', 6]]]

        (5+3)*6
        [[[5, '+', 3], '*', 6]]

        -2--11
        [[['-', 2], '-', ['-', 11]]]
    \"\"\"
    # captive version of FollowedBy that does not do parse actions or capture results names
    class _FB(FollowedBy):
        def parseImpl(self, instring, loc, doActions=True):
            self.expr.try_parse(instring, loc)
            return loc, []

    _FB.__name__ = \"FollowedBy>\"

    ret = Forward()
    if isinstance(lpar, str):
        lpar = Suppress(lpar)
    if isinstance(rpar, str):
        rpar = Suppress(rpar)

    # if lpar and rpar are not suppressed, wrap in group
    if not (isinstance(rpar, Suppress) and isinstance(rpar, Suppress)):
        lastExpr = base_expr | Group(lpar + ret + rpar)
    else:
        lastExpr = base_expr | (lpar + ret + rpar)

    for i, operDef in enumerate(op_list):
        opExpr, arity, rightLeftAssoc, pa = (operDef + (None,))[:4]
        if isinstance(opExpr, str_type):
            opExpr = ParserElement._literalStringClass(opExpr)
        if arity == 3:
            if not isinstance(opExpr, (tuple, list)) or len(opExpr) != 2:
                raise ValueError(
                    \"if numterms=3, opExpr must be a tuple or list of two expressions\"
                )
            opExpr1, opExpr2 = opExpr
            term_name = \"{}{} term\".format(opExpr1, opExpr2)
        else:
            term_name = \"{} term\".format(opExpr)

        if not 1 <= arity <= 3:
            raise ValueError(\"operator must be unary (1), binary (2), or ternary (3)\")

        if rightLeftAssoc not in (OpAssoc.LEFT, OpAssoc.RIGHT):
            raise ValueError(\"operator must indicate right or left associativity\")

        thisExpr: Forward = Forward().set_name(term_name)
        if rightLeftAssoc is OpAssoc.LEFT:
            if arity == 1:
                matchExpr = _FB(lastExpr + opExpr) + Group(lastExpr + opExpr[1, ...])
            elif arity == 2:
                if opExpr is not None:
                    matchExpr = _FB(lastExpr + opExpr + lastExpr) + Group(
                        lastExpr + (opExpr + lastExpr)[1, ...]
                    )
                else:
                    matchExpr = _FB(lastExpr + lastExpr) + Group(lastExpr[2, ...])
            elif arity == 3:
                matchExpr = _FB(
                    lastExpr + opExpr1 + lastExpr + opExpr2 + lastExpr
                ) + Group(lastExpr + OneOrMore(opExpr1 + lastExpr + opExpr2 + lastExpr))
        elif rightLeftAssoc is OpAssoc.RIGHT:
            if arity == 1:
                # try to avoid LR with this extra test
                if not isinstance(opExpr, Opt):
                    opExpr = Opt(opExpr)
                matchExpr = _FB(opExpr.expr + thisExpr) + Group(opExpr + thisExpr)
            elif arity == 2:
                if opExpr is not None:
                    matchExpr = _FB(lastExpr + opExpr + thisExpr) + Group(
                        lastExpr + (opExpr + thisExpr)[1, ...]
                    )
                else:
                    matchExpr = _FB(lastExpr + thisExpr) + Group(
                        lastExpr + thisExpr[1, ...]
                    )
            elif arity == 3:
                matchExpr = _FB(
                    lastExpr + opExpr1 + thisExpr + opExpr2 + thisExpr
                ) + Group(lastExpr + opExpr1 + thisExpr + opExpr2 + thisExpr)
        if pa:
            if isinstance(pa, (tuple, list)):
                matchExpr.set_parse_action(*pa)
            else:
                matchExpr.set_parse_action(pa)
        thisExpr <<= (matchExpr | lastExpr).setName(term_name)
        lastExpr = thisExpr
    ret <<= lastExpr
    return ret


def indentedBlock(blockStatementExpr, indentStack, indent=True, backup_stacks=[]):
    \"\"\"
    (DEPRECATED - use IndentedBlock class instead)
    Helper method for defining space-delimited indentation blocks,
    such as those used to define block statements in Python source code.

    Parameters:

    - ``blockStatementExpr`` - expression defining syntax of statement that
      is repeated within the indented block
    - ``indentStack`` - list created by caller to manage indentation stack
      (multiple ``statementWithIndentedBlock`` expressions within a single
      grammar should share a common ``indentStack``)
    - ``indent`` - boolean indicating whether block must be indented beyond
      the current level; set to ``False`` for block of left-most statements
      (default= ``True``)

    A valid block must contain at least one ``blockStatement``.

    (Note that indentedBlock uses internal parse actions which make it
    incompatible with packrat parsing.)

    Example::

        data = '''
        def A(z):
          A1
          B = 100
          G = A2
          A2
          A3
        B
        def BB(a,b,c):
          BB1
          def BBA():
            bba1
            bba2
            bba3
        C
        D
        def spam(x,y):
             def eggs(z):
                 pass
        '''


        indentStack = [1]
        stmt = Forward()

        identifier = Word(alphas, alphanums)
        funcDecl = (\"def\" + identifier + Group(\"(\" + Opt(delimitedList(identifier)) + \")\") + \":\")
        func_body = indentedBlock(stmt, indentStack)
        funcDef = Group(funcDecl + func_body)

        rvalue = Forward()
        funcCall = Group(identifier + \"(\" + Opt(delimitedList(rvalue)) + \")\")
        rvalue << (funcCall | identifier | Word(nums))
        assignment = Group(identifier + \"=\" + rvalue)
        stmt << (funcDef | assignment | identifier)

        module_body = stmt[1, ...]

        parseTree = module_body.parseString(data)
        parseTree.pprint()

    prints::

        [['def',
          'A',
          ['(', 'z', ')'],
          ':',
          [['A1'], [['B', '=', '100']], [['G', '=', 'A2']], ['A2'], ['A3']]],
         'B',
         ['def',
          'BB',
          ['(', 'a', 'b', 'c', ')'],
          ':',
          [['BB1'], [['def', 'BBA', ['(', ')'], ':', [['bba1'], ['bba2'], ['bba3']]]]]],
         'C',
         'D',
         ['def',
          'spam',
          ['(', 'x', 'y', ')'],
          ':',
          [[['def', 'eggs', ['(', 'z', ')'], ':', [['pass']]]]]]]
    \"\"\"
    backup_stacks.append(indentStack[:])

    def reset_stack():
        indentStack[:] = backup_stacks[-1]

    def checkPeerIndent(s, l, t):
        if l >= len(s):
            return
        curCol = col(l, s)
        if curCol != indentStack[-1]:
            if curCol > indentStack[-1]:
                raise ParseException(s, l, \"illegal nesting\")
            raise ParseException(s, l, \"not a peer entry\")

    def checkSubIndent(s, l, t):
        curCol = col(l, s)
        if curCol > indentStack[-1]:
            indentStack.append(curCol)
        else:
            raise ParseException(s, l, \"not a subentry\")

    def checkUnindent(s, l, t):
        if l >= len(s):
            return
        curCol = col(l, s)
        if not (indentStack and curCol in indentStack):
            raise ParseException(s, l, \"not an unindent\")
        if curCol < indentStack[-1]:
            indentStack.pop()

    NL = OneOrMore(LineEnd().set_whitespace_chars(\"\\t \").suppress())
    INDENT = (Empty() + Empty().set_parse_action(checkSubIndent)).set_name(\"INDENT\")
    PEER = Empty().set_parse_action(checkPeerIndent).set_name(\"\")
    UNDENT = Empty().set_parse_action(checkUnindent).set_name(\"UNINDENT\")
    if indent:
        smExpr = Group(
            Opt(NL)
            + INDENT
            + OneOrMore(PEER + Group(blockStatementExpr) + Opt(NL))
            + UNDENT
        )
    else:
        smExpr = Group(
            Opt(NL)
            + OneOrMore(PEER + Group(blockStatementExpr) + Opt(NL))
            + Opt(UNDENT)
        )

    # add a parse action to remove backup_stack from list of backups
    smExpr.add_parse_action(
        lambda: backup_stacks.pop(-1) and None if backup_stacks else None
    )
    smExpr.set_fail_action(lambda a, b, c, d: reset_stack())
    blockStatementExpr.ignore(_bslash + LineEnd())
    return smExpr.set_name(\"indented block\")


# it's easy to get these comment structures wrong - they're very common, so may as well make them available
c_style_comment = Combine(Regex(r\"/\\*(?:[^*]|\\*(?!/))*\") + \"*/\").set_name(
    \"C style comment\"
)
\"Comment of the form ``/* ... */``\"

html_comment = Regex(r\"<!--[\\s\\S]*?-->\").set_name(\"HTML comment\")
\"Comment of the form ``<!-- ... -->``\"

rest_of_line = Regex(r\".*\").leave_whitespace().set_name(\"rest of line\")
dbl_slash_comment = Regex(r\"//(?:\\\\\\n|[^\\n])*\").set_name(\"// comment\")
\"Comment of the form ``// ... (to end of line)``\"

cpp_style_comment = Combine(
    Regex(r\"/\\*(?:[^*]|\\*(?!/))*\") + \"*/\" | dbl_slash_comment
).set_name(\"C++ style comment\")
\"Comment of either form :class:`c_style_comment` or :class:`dbl_slash_comment`\"

java_style_comment = cpp_style_comment
\"Same as :class:`cpp_style_comment`\"

python_style_comment = Regex(r\"#.*\").set_name(\"Python style comment\")
\"Comment of the form ``# ... (to end of line)``\"


# build list of built-in expressions, for future reference if a global default value
# gets updated
_builtin_exprs: List[ParserElement] = [
    v for v in vars().values() if isinstance(v, ParserElement)
]


# pre-PEP8 compatible names
delimitedList = delimited_list
countedArray = counted_array
matchPreviousLiteral = match_previous_literal
matchPreviousExpr = match_previous_expr
oneOf = one_of
dictOf = dict_of
originalTextFor = original_text_for
nestedExpr = nested_expr
makeHTMLTags = make_html_tags
makeXMLTags = make_xml_tags
anyOpenTag, anyCloseTag = any_open_tag, any_close_tag
commonHTMLEntity = common_html_entity
replaceHTMLEntity = replace_html_entity
opAssoc = OpAssoc
infixNotation = infix_notation
cStyleComment = c_style_comment
htmlComment = html_comment
restOfLine = rest_of_line
dblSlashComment = dbl_slash_comment
cppStyleComment = cpp_style_comment
javaStyleComment = java_style_comment
pythonStyleComment = python_style_comment

"""
module_dict["pyparsing/testing.py"]="""
# testing.py

from contextlib import contextmanager
import typing

from .core import (
    ParserElement,
    ParseException,
    Keyword,
    __diag__,
    __compat__,
)


class pyparsing_test:
    \"\"\"
    namespace class for classes useful in writing unit tests
    \"\"\"

    class reset_pyparsing_context:
        \"\"\"
        Context manager to be used when writing unit tests that modify pyparsing config values:
        - packrat parsing
        - bounded recursion parsing
        - default whitespace characters.
        - default keyword characters
        - literal string auto-conversion class
        - __diag__ settings

        Example::

            with reset_pyparsing_context():
                # test that literals used to construct a grammar are automatically suppressed
                ParserElement.inlineLiteralsUsing(Suppress)

                term = Word(alphas) | Word(nums)
                group = Group('(' + term[...] + ')')

                # assert that the '()' characters are not included in the parsed tokens
                self.assertParseAndCheckList(group, \"(abc 123 def)\", ['abc', '123', 'def'])

            # after exiting context manager, literals are converted to Literal expressions again
        \"\"\"

        def __init__(self):
            self._save_context = {}

        def save(self):
            self._save_context[\"default_whitespace\"] = ParserElement.DEFAULT_WHITE_CHARS
            self._save_context[\"default_keyword_chars\"] = Keyword.DEFAULT_KEYWORD_CHARS

            self._save_context[
                \"literal_string_class\"
            ] = ParserElement._literalStringClass

            self._save_context[\"verbose_stacktrace\"] = ParserElement.verbose_stacktrace

            self._save_context[\"packrat_enabled\"] = ParserElement._packratEnabled
            if ParserElement._packratEnabled:
                self._save_context[
                    \"packrat_cache_size\"
                ] = ParserElement.packrat_cache.size
            else:
                self._save_context[\"packrat_cache_size\"] = None
            self._save_context[\"packrat_parse\"] = ParserElement._parse
            self._save_context[
                \"recursion_enabled\"
            ] = ParserElement._left_recursion_enabled

            self._save_context[\"__diag__\"] = {
                name: getattr(__diag__, name) for name in __diag__._all_names
            }

            self._save_context[\"__compat__\"] = {
                \"collect_all_And_tokens\": __compat__.collect_all_And_tokens
            }

            return self

        def restore(self):
            # reset pyparsing global state
            if (
                ParserElement.DEFAULT_WHITE_CHARS
                != self._save_context[\"default_whitespace\"]
            ):
                ParserElement.set_default_whitespace_chars(
                    self._save_context[\"default_whitespace\"]
                )

            ParserElement.verbose_stacktrace = self._save_context[\"verbose_stacktrace\"]

            Keyword.DEFAULT_KEYWORD_CHARS = self._save_context[\"default_keyword_chars\"]
            ParserElement.inlineLiteralsUsing(
                self._save_context[\"literal_string_class\"]
            )

            for name, value in self._save_context[\"__diag__\"].items():
                (__diag__.enable if value else __diag__.disable)(name)

            ParserElement._packratEnabled = False
            if self._save_context[\"packrat_enabled\"]:
                ParserElement.enable_packrat(self._save_context[\"packrat_cache_size\"])
            else:
                ParserElement._parse = self._save_context[\"packrat_parse\"]
            ParserElement._left_recursion_enabled = self._save_context[
                \"recursion_enabled\"
            ]

            __compat__.collect_all_And_tokens = self._save_context[\"__compat__\"]

            return self

        def copy(self):
            ret = type(self)()
            ret._save_context.update(self._save_context)
            return ret

        def __enter__(self):
            return self.save()

        def __exit__(self, *args):
            self.restore()

    class TestParseResultsAsserts:
        \"\"\"
        A mixin class to add parse results assertion methods to normal unittest.TestCase classes.
        \"\"\"

        def assertParseResultsEquals(
            self, result, expected_list=None, expected_dict=None, msg=None
        ):
            \"\"\"
            Unit test assertion to compare a :class:`ParseResults` object with an optional ``expected_list``,
            and compare any defined results names with an optional ``expected_dict``.
            \"\"\"
            if expected_list is not None:
                self.assertEqual(expected_list, result.as_list(), msg=msg)
            if expected_dict is not None:
                self.assertEqual(expected_dict, result.as_dict(), msg=msg)

        def assertParseAndCheckList(
            self, expr, test_string, expected_list, msg=None, verbose=True
        ):
            \"\"\"
            Convenience wrapper assert to test a parser element and input string, and assert that
            the resulting ``ParseResults.asList()`` is equal to the ``expected_list``.
            \"\"\"
            result = expr.parse_string(test_string, parse_all=True)
            if verbose:
                print(result.dump())
            else:
                print(result.as_list())
            self.assertParseResultsEquals(result, expected_list=expected_list, msg=msg)

        def assertParseAndCheckDict(
            self, expr, test_string, expected_dict, msg=None, verbose=True
        ):
            \"\"\"
            Convenience wrapper assert to test a parser element and input string, and assert that
            the resulting ``ParseResults.asDict()`` is equal to the ``expected_dict``.
            \"\"\"
            result = expr.parse_string(test_string, parseAll=True)
            if verbose:
                print(result.dump())
            else:
                print(result.as_list())
            self.assertParseResultsEquals(result, expected_dict=expected_dict, msg=msg)

        def assertRunTestResults(
            self, run_tests_report, expected_parse_results=None, msg=None
        ):
            \"\"\"
            Unit test assertion to evaluate output of ``ParserElement.runTests()``. If a list of
            list-dict tuples is given as the ``expected_parse_results`` argument, then these are zipped
            with the report tuples returned by ``runTests`` and evaluated using ``assertParseResultsEquals``.
            Finally, asserts that the overall ``runTests()`` success value is ``True``.

            :param run_tests_report: tuple(bool, [tuple(str, ParseResults or Exception)]) returned from runTests
            :param expected_parse_results (optional): [tuple(str, list, dict, Exception)]
            \"\"\"
            run_test_success, run_test_results = run_tests_report

            if expected_parse_results is not None:
                merged = [
                    (*rpt, expected)
                    for rpt, expected in zip(run_test_results, expected_parse_results)
                ]
                for test_string, result, expected in merged:
                    # expected should be a tuple containing a list and/or a dict or an exception,
                    # and optional failure message string
                    # an empty tuple will skip any result validation
                    fail_msg = next(
                        (exp for exp in expected if isinstance(exp, str)), None
                    )
                    expected_exception = next(
                        (
                            exp
                            for exp in expected
                            if isinstance(exp, type) and issubclass(exp, Exception)
                        ),
                        None,
                    )
                    if expected_exception is not None:
                        with self.assertRaises(
                            expected_exception=expected_exception, msg=fail_msg or msg
                        ):
                            if isinstance(result, Exception):
                                raise result
                    else:
                        expected_list = next(
                            (exp for exp in expected if isinstance(exp, list)), None
                        )
                        expected_dict = next(
                            (exp for exp in expected if isinstance(exp, dict)), None
                        )
                        if (expected_list, expected_dict) != (None, None):
                            self.assertParseResultsEquals(
                                result,
                                expected_list=expected_list,
                                expected_dict=expected_dict,
                                msg=fail_msg or msg,
                            )
                        else:
                            # warning here maybe?
                            print(\"no validation for {!r}\".format(test_string))

            # do this last, in case some specific test results can be reported instead
            self.assertTrue(
                run_test_success, msg=msg if msg is not None else \"failed runTests\"
            )

        @contextmanager
        def assertRaisesParseException(self, exc_type=ParseException, msg=None):
            with self.assertRaises(exc_type, msg=msg):
                yield

    @staticmethod
    def with_line_numbers(
        s: str,
        start_line: typing.Optional[int] = None,
        end_line: typing.Optional[int] = None,
        expand_tabs: bool = True,
        eol_mark: str = \"|\",
        mark_spaces: typing.Optional[str] = None,
        mark_control: typing.Optional[str] = None,
    ) -> str:
        \"\"\"
        Helpful method for debugging a parser - prints a string with line and column numbers.
        (Line and column numbers are 1-based.)

        :param s: tuple(bool, str - string to be printed with line and column numbers
        :param start_line: int - (optional) starting line number in s to print (default=1)
        :param end_line: int - (optional) ending line number in s to print (default=len(s))
        :param expand_tabs: bool - (optional) expand tabs to spaces, to match the pyparsing default
        :param eol_mark: str - (optional) string to mark the end of lines, helps visualize trailing spaces (default=\"|\")
        :param mark_spaces: str - (optional) special character to display in place of spaces
        :param mark_control: str - (optional) convert non-printing control characters to a placeholding
                                 character; valid values:
                                 - \"unicode\" - replaces control chars with Unicode symbols, such as \"␍\" and \"␊\"
                                 - any single character string - replace control characters with given string
                                 - None (default) - string is displayed as-is

        :return: str - input string with leading line numbers and column number headers
        \"\"\"
        if expand_tabs:
            s = s.expandtabs()
        if mark_control is not None:
            if mark_control == \"unicode\":
                tbl = str.maketrans(
                    {c: u for c, u in zip(range(0, 33), range(0x2400, 0x2433))}
                    | {127: 0x2421}
                )
                eol_mark = \"\"
            else:
                tbl = str.maketrans(
                    {c: mark_control for c in list(range(0, 32)) + [127]}
                )
            s = s.translate(tbl)
        if mark_spaces is not None and mark_spaces != \" \":
            if mark_spaces == \"unicode\":
                tbl = str.maketrans({9: 0x2409, 32: 0x2423})
                s = s.translate(tbl)
            else:
                s = s.replace(\" \", mark_spaces)
        if start_line is None:
            start_line = 1
        if end_line is None:
            end_line = len(s)
        end_line = min(end_line, len(s))
        start_line = min(max(1, start_line), end_line)

        if mark_control != \"unicode\":
            s_lines = s.splitlines()[start_line - 1 : end_line]
        else:
            s_lines = [line + \"␊\" for line in s.split(\"␊\")[start_line - 1 : end_line]]
        if not s_lines:
            return \"\"

        lineno_width = len(str(end_line))
        max_line_len = max(len(line) for line in s_lines)
        lead = \" \" * (lineno_width + 1)
        if max_line_len >= 99:
            header0 = (
                lead
                + \"\".join(
                    \"{}{}\".format(\" \" * 99, (i + 1) % 100)
                    for i in range(max(max_line_len // 100, 1))
                )
                + \"\\n\"
            )
        else:
            header0 = \"\"
        header1 = (
            header0
            + lead
            + \"\".join(
                \"         {}\".format((i + 1) % 10)
                for i in range(-(-max_line_len // 10))
            )
            + \"\\n\"
        )
        header2 = lead + \"1234567890\" * (-(-max_line_len // 10)) + \"\\n\"
        return (
            header1
            + header2
            + \"\\n\".join(
                \"{:{}d}:{}{}\".format(i, lineno_width, line, eol_mark)
                for i, line in enumerate(s_lines, start=start_line)
            )
            + \"\\n\"
        )

"""
module_dict["pyparsing/exceptions.py"]="""
# exceptions.py

import re
import sys
import typing

from .util import col, line, lineno, _collapse_string_to_ranges
from .unicode import pyparsing_unicode as ppu


class ExceptionWordUnicode(ppu.Latin1, ppu.LatinA, ppu.LatinB, ppu.Greek, ppu.Cyrillic):
    pass


_extract_alphanums = _collapse_string_to_ranges(ExceptionWordUnicode.alphanums)
_exception_word_extractor = re.compile(\"([\" + _extract_alphanums + \"]{1,16})|.\")


class ParseBaseException(Exception):
    \"\"\"base exception class for all parsing runtime exceptions\"\"\"

    # Performance tuning: we construct a *lot* of these, so keep this
    # constructor as small and fast as possible
    def __init__(
        self,
        pstr: str,
        loc: int = 0,
        msg: typing.Optional[str] = None,
        elem=None,
    ):
        self.loc = loc
        if msg is None:
            self.msg = pstr
            self.pstr = \"\"
        else:
            self.msg = msg
            self.pstr = pstr
        self.parser_element = self.parserElement = elem
        self.args = (pstr, loc, msg)

    @staticmethod
    def explain_exception(exc, depth=16):
        \"\"\"
        Method to take an exception and translate the Python internal traceback into a list
        of the pyparsing expressions that caused the exception to be raised.

        Parameters:

        - exc - exception raised during parsing (need not be a ParseException, in support
          of Python exceptions that might be raised in a parse action)
        - depth (default=16) - number of levels back in the stack trace to list expression
          and function names; if None, the full stack trace names will be listed; if 0, only
          the failing input line, marker, and exception string will be shown

        Returns a multi-line string listing the ParserElements and/or function names in the
        exception's stack trace.
        \"\"\"
        import inspect
        from .core import ParserElement

        if depth is None:
            depth = sys.getrecursionlimit()
        ret = []
        if isinstance(exc, ParseBaseException):
            ret.append(exc.line)
            ret.append(\" \" * (exc.column - 1) + \"^\")
        ret.append(\"{}: {}\".format(type(exc).__name__, exc))

        if depth > 0:
            callers = inspect.getinnerframes(exc.__traceback__, context=depth)
            seen = set()
            for i, ff in enumerate(callers[-depth:]):
                frm = ff[0]

                f_self = frm.f_locals.get(\"self\", None)
                if isinstance(f_self, ParserElement):
                    if frm.f_code.co_name not in (\"parseImpl\", \"_parseNoCache\"):
                        continue
                    if id(f_self) in seen:
                        continue
                    seen.add(id(f_self))

                    self_type = type(f_self)
                    ret.append(
                        \"{}.{} - {}\".format(
                            self_type.__module__, self_type.__name__, f_self
                        )
                    )

                elif f_self is not None:
                    self_type = type(f_self)
                    ret.append(\"{}.{}\".format(self_type.__module__, self_type.__name__))

                else:
                    code = frm.f_code
                    if code.co_name in (\"wrapper\", \"<module>\"):
                        continue

                    ret.append(\"{}\".format(code.co_name))

                depth -= 1
                if not depth:
                    break

        return \"\\n\".join(ret)

    @classmethod
    def _from_exception(cls, pe):
        \"\"\"
        internal factory method to simplify creating one type of ParseException
        from another - avoids having __init__ signature conflicts among subclasses
        \"\"\"
        return cls(pe.pstr, pe.loc, pe.msg, pe.parserElement)

    @property
    def line(self) -> str:
        \"\"\"
        Return the line of text where the exception occurred.
        \"\"\"
        return line(self.loc, self.pstr)

    @property
    def lineno(self) -> int:
        \"\"\"
        Return the 1-based line number of text where the exception occurred.
        \"\"\"
        return lineno(self.loc, self.pstr)

    @property
    def col(self) -> int:
        \"\"\"
        Return the 1-based column on the line of text where the exception occurred.
        \"\"\"
        return col(self.loc, self.pstr)

    @property
    def column(self) -> int:
        \"\"\"
        Return the 1-based column on the line of text where the exception occurred.
        \"\"\"
        return col(self.loc, self.pstr)

    def __str__(self) -> str:
        if self.pstr:
            if self.loc >= len(self.pstr):
                foundstr = \", found end of text\"
            else:
                # pull out next word at error location
                found_match = _exception_word_extractor.match(self.pstr, self.loc)
                if found_match is not None:
                    found = found_match.group(0)
                else:
                    found = self.pstr[self.loc : self.loc + 1]
                foundstr = (\", found %r\" % found).replace(r\"\\\\\", \"\\\\\")
        else:
            foundstr = \"\"
        return \"{}{}  (at char {}), (line:{}, col:{})\".format(
            self.msg, foundstr, self.loc, self.lineno, self.column
        )

    def __repr__(self):
        return str(self)

    def mark_input_line(self, marker_string: str = None, *, markerString=\">!<\") -> str:
        \"\"\"
        Extracts the exception line from the input string, and marks
        the location of the exception with a special symbol.
        \"\"\"
        markerString = marker_string if marker_string is not None else markerString
        line_str = self.line
        line_column = self.column - 1
        if markerString:
            line_str = \"\".join(
                (line_str[:line_column], markerString, line_str[line_column:])
            )
        return line_str.strip()

    def explain(self, depth=16) -> str:
        \"\"\"
        Method to translate the Python internal traceback into a list
        of the pyparsing expressions that caused the exception to be raised.

        Parameters:

        - depth (default=16) - number of levels back in the stack trace to list expression
          and function names; if None, the full stack trace names will be listed; if 0, only
          the failing input line, marker, and exception string will be shown

        Returns a multi-line string listing the ParserElements and/or function names in the
        exception's stack trace.

        Example::

            expr = pp.Word(pp.nums) * 3
            try:
                expr.parse_string(\"123 456 A789\")
            except pp.ParseException as pe:
                print(pe.explain(depth=0))

        prints::

            123 456 A789
                    ^
            ParseException: Expected W:(0-9), found 'A'  (at char 8), (line:1, col:9)

        Note: the diagnostic output will include string representations of the expressions
        that failed to parse. These representations will be more helpful if you use `set_name` to
        give identifiable names to your expressions. Otherwise they will use the default string
        forms, which may be cryptic to read.

        Note: pyparsing's default truncation of exception tracebacks may also truncate the
        stack of expressions that are displayed in the ``explain`` output. To get the full listing
        of parser expressions, you may have to set ``ParserElement.verbose_stacktrace = True``
        \"\"\"
        return self.explain_exception(self, depth)

    markInputline = mark_input_line


class ParseException(ParseBaseException):
    \"\"\"
    Exception thrown when a parse expression doesn't match the input string

    Example::

        try:
            Word(nums).set_name(\"integer\").parse_string(\"ABC\")
        except ParseException as pe:
            print(pe)
            print(\"column: {}\".format(pe.column))

    prints::

       Expected integer (at char 0), (line:1, col:1)
        column: 1

    \"\"\"


class ParseFatalException(ParseBaseException):
    \"\"\"
    User-throwable exception thrown when inconsistent parse content
    is found; stops all parsing immediately
    \"\"\"


class ParseSyntaxException(ParseFatalException):
    \"\"\"
    Just like :class:`ParseFatalException`, but thrown internally
    when an :class:`ErrorStop<And._ErrorStop>` ('-' operator) indicates
    that parsing is to stop immediately because an unbacktrackable
    syntax error has been found.
    \"\"\"


class RecursiveGrammarException(Exception):
    \"\"\"
    Exception thrown by :class:`ParserElement.validate` if the
    grammar could be left-recursive; parser may need to enable
    left recursion using :class:`ParserElement.enable_left_recursion<ParserElement.enable_left_recursion>`
    \"\"\"

    def __init__(self, parseElementList):
        self.parseElementTrace = parseElementList

    def __str__(self) -> str:
        return \"RecursiveGrammarException: {}\".format(self.parseElementTrace)

"""
module_dict["pyparsing/results.py"]="""
# results.py
from collections.abc import MutableMapping, Mapping, MutableSequence, Iterator
import pprint
from weakref import ref as wkref
from typing import Tuple, Any

str_type: Tuple[type, ...] = (str, bytes)
_generator_type = type((_ for _ in ()))


class _ParseResultsWithOffset:
    __slots__ = [\"tup\"]

    def __init__(self, p1, p2):
        self.tup = (p1, p2)

    def __getitem__(self, i):
        return self.tup[i]

    def __getstate__(self):
        return self.tup

    def __setstate__(self, *args):
        self.tup = args[0]


class ParseResults:
    \"\"\"Structured parse results, to provide multiple means of access to
    the parsed data:

    - as a list (``len(results)``)
    - by list index (``results[0], results[1]``, etc.)
    - by attribute (``results.<results_name>`` - see :class:`ParserElement.set_results_name`)

    Example::

        integer = Word(nums)
        date_str = (integer.set_results_name(\"year\") + '/'
                    + integer.set_results_name(\"month\") + '/'
                    + integer.set_results_name(\"day\"))
        # equivalent form:
        # date_str = (integer(\"year\") + '/'
        #             + integer(\"month\") + '/'
        #             + integer(\"day\"))

        # parse_string returns a ParseResults object
        result = date_str.parse_string(\"1999/12/31\")

        def test(s, fn=repr):
            print(\"{} -> {}\".format(s, fn(eval(s))))
        test(\"list(result)\")
        test(\"result[0]\")
        test(\"result['month']\")
        test(\"result.day\")
        test(\"'month' in result\")
        test(\"'minutes' in result\")
        test(\"result.dump()\", str)

    prints::

        list(result) -> ['1999', '/', '12', '/', '31']
        result[0] -> '1999'
        result['month'] -> '12'
        result.day -> '31'
        'month' in result -> True
        'minutes' in result -> False
        result.dump() -> ['1999', '/', '12', '/', '31']
        - day: '31'
        - month: '12'
        - year: '1999'
    \"\"\"

    _null_values: Tuple[Any, ...] = (None, [], \"\", ())

    __slots__ = [
        \"_name\",
        \"_parent\",
        \"_all_names\",
        \"_modal\",
        \"_toklist\",
        \"_tokdict\",
        \"__weakref__\",
    ]

    class List(list):
        \"\"\"
        Simple wrapper class to distinguish parsed list results that should be preserved
        as actual Python lists, instead of being converted to :class:`ParseResults`:

            LBRACK, RBRACK = map(pp.Suppress, \"[]\")
            element = pp.Forward()
            item = ppc.integer
            element_list = LBRACK + pp.delimited_list(element) + RBRACK

            # add parse actions to convert from ParseResults to actual Python collection types
            def as_python_list(t):
                return pp.ParseResults.List(t.as_list())
            element_list.add_parse_action(as_python_list)

            element <<= item | element_list

            element.run_tests('''
                100
                [2,3,4]
                [[2, 1],3,4]
                [(2, 1),3,4]
                (2,3,4)
                ''', post_parse=lambda s, r: (r[0], type(r[0])))

        prints:

            100
            (100, <class 'int'>)

            [2,3,4]
            ([2, 3, 4], <class 'list'>)

            [[2, 1],3,4]
            ([[2, 1], 3, 4], <class 'list'>)

        (Used internally by :class:`Group` when `aslist=True`.)
        \"\"\"

        def __new__(cls, contained=None):
            if contained is None:
                contained = []

            if not isinstance(contained, list):
                raise TypeError(
                    \"{} may only be constructed with a list,\"
                    \" not {}\".format(cls.__name__, type(contained).__name__)
                )

            return list.__new__(cls)

    def __new__(cls, toklist=None, name=None, **kwargs):
        if isinstance(toklist, ParseResults):
            return toklist
        self = object.__new__(cls)
        self._name = None
        self._parent = None
        self._all_names = set()

        if toklist is None:
            self._toklist = []
        elif isinstance(toklist, (list, _generator_type)):
            self._toklist = (
                [toklist[:]]
                if isinstance(toklist, ParseResults.List)
                else list(toklist)
            )
        else:
            self._toklist = [toklist]
        self._tokdict = dict()
        return self

    # Performance tuning: we construct a *lot* of these, so keep this
    # constructor as small and fast as possible
    def __init__(
        self, toklist=None, name=None, asList=True, modal=True, isinstance=isinstance
    ):
        self._modal = modal
        if name is not None and name != \"\":
            if isinstance(name, int):
                name = str(name)
            if not modal:
                self._all_names = {name}
            self._name = name
            if toklist not in self._null_values:
                if isinstance(toklist, (str_type, type)):
                    toklist = [toklist]
                if asList:
                    if isinstance(toklist, ParseResults):
                        self[name] = _ParseResultsWithOffset(
                            ParseResults(toklist._toklist), 0
                        )
                    else:
                        self[name] = _ParseResultsWithOffset(
                            ParseResults(toklist[0]), 0
                        )
                    self[name]._name = name
                else:
                    try:
                        self[name] = toklist[0]
                    except (KeyError, TypeError, IndexError):
                        if toklist is not self:
                            self[name] = toklist
                        else:
                            self._name = name

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self._toklist[i]
        else:
            if i not in self._all_names:
                return self._tokdict[i][-1][0]
            else:
                return ParseResults([v[0] for v in self._tokdict[i]])

    def __setitem__(self, k, v, isinstance=isinstance):
        if isinstance(v, _ParseResultsWithOffset):
            self._tokdict[k] = self._tokdict.get(k, list()) + [v]
            sub = v[0]
        elif isinstance(k, (int, slice)):
            self._toklist[k] = v
            sub = v
        else:
            self._tokdict[k] = self._tokdict.get(k, list()) + [
                _ParseResultsWithOffset(v, 0)
            ]
            sub = v
        if isinstance(sub, ParseResults):
            sub._parent = wkref(self)

    def __delitem__(self, i):
        if isinstance(i, (int, slice)):
            mylen = len(self._toklist)
            del self._toklist[i]

            # convert int to slice
            if isinstance(i, int):
                if i < 0:
                    i += mylen
                i = slice(i, i + 1)
            # get removed indices
            removed = list(range(*i.indices(mylen)))
            removed.reverse()
            # fixup indices in token dictionary
            for name, occurrences in self._tokdict.items():
                for j in removed:
                    for k, (value, position) in enumerate(occurrences):
                        occurrences[k] = _ParseResultsWithOffset(
                            value, position - (position > j)
                        )
        else:
            del self._tokdict[i]

    def __contains__(self, k) -> bool:
        return k in self._tokdict

    def __len__(self) -> int:
        return len(self._toklist)

    def __bool__(self) -> bool:
        return not not (self._toklist or self._tokdict)

    def __iter__(self) -> Iterator:
        return iter(self._toklist)

    def __reversed__(self) -> Iterator:
        return iter(self._toklist[::-1])

    def keys(self):
        return iter(self._tokdict)

    def values(self):
        return (self[k] for k in self.keys())

    def items(self):
        return ((k, self[k]) for k in self.keys())

    def haskeys(self) -> bool:
        \"\"\"
        Since ``keys()`` returns an iterator, this method is helpful in bypassing
        code that looks for the existence of any defined results names.\"\"\"
        return bool(self._tokdict)

    def pop(self, *args, **kwargs):
        \"\"\"
        Removes and returns item at specified index (default= ``last``).
        Supports both ``list`` and ``dict`` semantics for ``pop()``. If
        passed no argument or an integer argument, it will use ``list``
        semantics and pop tokens from the list of parsed tokens. If passed
        a non-integer argument (most likely a string), it will use ``dict``
        semantics and pop the corresponding value from any defined results
        names. A second default return value argument is supported, just as in
        ``dict.pop()``.

        Example::

            numlist = Word(nums)[...]
            print(numlist.parse_string(\"0 123 321\")) # -> ['0', '123', '321']

            def remove_first(tokens):
                tokens.pop(0)
            numlist.add_parse_action(remove_first)
            print(numlist.parse_string(\"0 123 321\")) # -> ['123', '321']

            label = Word(alphas)
            patt = label(\"LABEL\") + Word(nums)[1, ...]
            print(patt.parse_string(\"AAB 123 321\").dump())

            # Use pop() in a parse action to remove named result (note that corresponding value is not
            # removed from list form of results)
            def remove_LABEL(tokens):
                tokens.pop(\"LABEL\")
                return tokens
            patt.add_parse_action(remove_LABEL)
            print(patt.parse_string(\"AAB 123 321\").dump())

        prints::

            ['AAB', '123', '321']
            - LABEL: 'AAB'

            ['AAB', '123', '321']
        \"\"\"
        if not args:
            args = [-1]
        for k, v in kwargs.items():
            if k == \"default\":
                args = (args[0], v)
            else:
                raise TypeError(
                    \"pop() got an unexpected keyword argument {!r}\".format(k)
                )
        if isinstance(args[0], int) or len(args) == 1 or args[0] in self:
            index = args[0]
            ret = self[index]
            del self[index]
            return ret
        else:
            defaultvalue = args[1]
            return defaultvalue

    def get(self, key, default_value=None):
        \"\"\"
        Returns named result matching the given key, or if there is no
        such name, then returns the given ``default_value`` or ``None`` if no
        ``default_value`` is specified.

        Similar to ``dict.get()``.

        Example::

            integer = Word(nums)
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")

            result = date_str.parse_string(\"1999/12/31\")
            print(result.get(\"year\")) # -> '1999'
            print(result.get(\"hour\", \"not specified\")) # -> 'not specified'
            print(result.get(\"hour\")) # -> None
        \"\"\"
        if key in self:
            return self[key]
        else:
            return default_value

    def insert(self, index, ins_string):
        \"\"\"
        Inserts new element at location index in the list of parsed tokens.

        Similar to ``list.insert()``.

        Example::

            numlist = Word(nums)[...]
            print(numlist.parse_string(\"0 123 321\")) # -> ['0', '123', '321']

            # use a parse action to insert the parse location in the front of the parsed results
            def insert_locn(locn, tokens):
                tokens.insert(0, locn)
            numlist.add_parse_action(insert_locn)
            print(numlist.parse_string(\"0 123 321\")) # -> [0, '0', '123', '321']
        \"\"\"
        self._toklist.insert(index, ins_string)
        # fixup indices in token dictionary
        for name, occurrences in self._tokdict.items():
            for k, (value, position) in enumerate(occurrences):
                occurrences[k] = _ParseResultsWithOffset(
                    value, position + (position > index)
                )

    def append(self, item):
        \"\"\"
        Add single element to end of ``ParseResults`` list of elements.

        Example::

            numlist = Word(nums)[...]
            print(numlist.parse_string(\"0 123 321\")) # -> ['0', '123', '321']

            # use a parse action to compute the sum of the parsed integers, and add it to the end
            def append_sum(tokens):
                tokens.append(sum(map(int, tokens)))
            numlist.add_parse_action(append_sum)
            print(numlist.parse_string(\"0 123 321\")) # -> ['0', '123', '321', 444]
        \"\"\"
        self._toklist.append(item)

    def extend(self, itemseq):
        \"\"\"
        Add sequence of elements to end of ``ParseResults`` list of elements.

        Example::

            patt = Word(alphas)[1, ...]

            # use a parse action to append the reverse of the matched strings, to make a palindrome
            def make_palindrome(tokens):
                tokens.extend(reversed([t[::-1] for t in tokens]))
                return ''.join(tokens)
            patt.add_parse_action(make_palindrome)
            print(patt.parse_string(\"lskdj sdlkjf lksd\")) # -> 'lskdjsdlkjflksddsklfjkldsjdksl'
        \"\"\"
        if isinstance(itemseq, ParseResults):
            self.__iadd__(itemseq)
        else:
            self._toklist.extend(itemseq)

    def clear(self):
        \"\"\"
        Clear all elements and results names.
        \"\"\"
        del self._toklist[:]
        self._tokdict.clear()

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            if name.startswith(\"__\"):
                raise AttributeError(name)
            return \"\"

    def __add__(self, other) -> \"ParseResults\":
        ret = self.copy()
        ret += other
        return ret

    def __iadd__(self, other) -> \"ParseResults\":
        if other._tokdict:
            offset = len(self._toklist)
            addoffset = lambda a: offset if a < 0 else a + offset
            otheritems = other._tokdict.items()
            otherdictitems = [
                (k, _ParseResultsWithOffset(v[0], addoffset(v[1])))
                for k, vlist in otheritems
                for v in vlist
            ]
            for k, v in otherdictitems:
                self[k] = v
                if isinstance(v[0], ParseResults):
                    v[0]._parent = wkref(self)

        self._toklist += other._toklist
        self._all_names |= other._all_names
        return self

    def __radd__(self, other) -> \"ParseResults\":
        if isinstance(other, int) and other == 0:
            # useful for merging many ParseResults using sum() builtin
            return self.copy()
        else:
            # this may raise a TypeError - so be it
            return other + self

    def __repr__(self) -> str:
        return \"{}({!r}, {})\".format(type(self).__name__, self._toklist, self.as_dict())

    def __str__(self) -> str:
        return (
            \"[\"
            + \", \".join(
                [
                    str(i) if isinstance(i, ParseResults) else repr(i)
                    for i in self._toklist
                ]
            )
            + \"]\"
        )

    def _asStringList(self, sep=\"\"):
        out = []
        for item in self._toklist:
            if out and sep:
                out.append(sep)
            if isinstance(item, ParseResults):
                out += item._asStringList()
            else:
                out.append(str(item))
        return out

    def as_list(self) -> list:
        \"\"\"
        Returns the parse results as a nested list of matching tokens, all converted to strings.

        Example::

            patt = Word(alphas)[1, ...]
            result = patt.parse_string(\"sldkj lsdkj sldkj\")
            # even though the result prints in string-like form, it is actually a pyparsing ParseResults
            print(type(result), result) # -> <class 'pyparsing.ParseResults'> ['sldkj', 'lsdkj', 'sldkj']

            # Use as_list() to create an actual list
            result_list = result.as_list()
            print(type(result_list), result_list) # -> <class 'list'> ['sldkj', 'lsdkj', 'sldkj']
        \"\"\"
        return [
            res.as_list() if isinstance(res, ParseResults) else res
            for res in self._toklist
        ]

    def as_dict(self) -> dict:
        \"\"\"
        Returns the named parse results as a nested dictionary.

        Example::

            integer = Word(nums)
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")

            result = date_str.parse_string('12/31/1999')
            print(type(result), repr(result)) # -> <class 'pyparsing.ParseResults'> (['12', '/', '31', '/', '1999'], {'day': [('1999', 4)], 'year': [('12', 0)], 'month': [('31', 2)]})

            result_dict = result.as_dict()
            print(type(result_dict), repr(result_dict)) # -> <class 'dict'> {'day': '1999', 'year': '12', 'month': '31'}

            # even though a ParseResults supports dict-like access, sometime you just need to have a dict
            import json
            print(json.dumps(result)) # -> Exception: TypeError: ... is not JSON serializable
            print(json.dumps(result.as_dict())) # -> {\"month\": \"31\", \"day\": \"1999\", \"year\": \"12\"}
        \"\"\"

        def to_item(obj):
            if isinstance(obj, ParseResults):
                return obj.as_dict() if obj.haskeys() else [to_item(v) for v in obj]
            else:
                return obj

        return dict((k, to_item(v)) for k, v in self.items())

    def copy(self) -> \"ParseResults\":
        \"\"\"
        Returns a new copy of a :class:`ParseResults` object.
        \"\"\"
        ret = ParseResults(self._toklist)
        ret._tokdict = self._tokdict.copy()
        ret._parent = self._parent
        ret._all_names |= self._all_names
        ret._name = self._name
        return ret

    def get_name(self):
        r\"\"\"
        Returns the results name for this token expression. Useful when several
        different expressions might match at a particular location.

        Example::

            integer = Word(nums)
            ssn_expr = Regex(r\"\\d\\d\\d-\\d\\d-\\d\\d\\d\\d\")
            house_number_expr = Suppress('#') + Word(nums, alphanums)
            user_data = (Group(house_number_expr)(\"house_number\")
                        | Group(ssn_expr)(\"ssn\")
                        | Group(integer)(\"age\"))
            user_info = user_data[1, ...]

            result = user_info.parse_string(\"22 111-22-3333 #221B\")
            for item in result:
                print(item.get_name(), ':', item[0])

        prints::

            age : 22
            ssn : 111-22-3333
            house_number : 221B
        \"\"\"
        if self._name:
            return self._name
        elif self._parent:
            par = self._parent()

            def find_in_parent(sub):
                return next(
                    (
                        k
                        for k, vlist in par._tokdict.items()
                        for v, loc in vlist
                        if sub is v
                    ),
                    None,
                )

            return find_in_parent(self) if par else None
        elif (
            len(self) == 1
            and len(self._tokdict) == 1
            and next(iter(self._tokdict.values()))[0][1] in (0, -1)
        ):
            return next(iter(self._tokdict.keys()))
        else:
            return None

    def dump(self, indent=\"\", full=True, include_list=True, _depth=0) -> str:
        \"\"\"
        Diagnostic method for listing out the contents of
        a :class:`ParseResults`. Accepts an optional ``indent`` argument so
        that this string can be embedded in a nested display of other data.

        Example::

            integer = Word(nums)
            date_str = integer(\"year\") + '/' + integer(\"month\") + '/' + integer(\"day\")

            result = date_str.parse_string('1999/12/31')
            print(result.dump())

        prints::

            ['1999', '/', '12', '/', '31']
            - day: '31'
            - month: '12'
            - year: '1999'
        \"\"\"
        out = []
        NL = \"\\n\"
        out.append(indent + str(self.as_list()) if include_list else \"\")

        if full:
            if self.haskeys():
                items = sorted((str(k), v) for k, v in self.items())
                for k, v in items:
                    if out:
                        out.append(NL)
                    out.append(\"{}{}- {}: \".format(indent, (\"  \" * _depth), k))
                    if isinstance(v, ParseResults):
                        if v:
                            out.append(
                                v.dump(
                                    indent=indent,
                                    full=full,
                                    include_list=include_list,
                                    _depth=_depth + 1,
                                )
                            )
                        else:
                            out.append(str(v))
                    else:
                        out.append(repr(v))
            if any(isinstance(vv, ParseResults) for vv in self):
                v = self
                for i, vv in enumerate(v):
                    if isinstance(vv, ParseResults):
                        out.append(
                            \"\\n{}{}[{}]:\\n{}{}{}\".format(
                                indent,
                                (\"  \" * (_depth)),
                                i,
                                indent,
                                (\"  \" * (_depth + 1)),
                                vv.dump(
                                    indent=indent,
                                    full=full,
                                    include_list=include_list,
                                    _depth=_depth + 1,
                                ),
                            )
                        )
                    else:
                        out.append(
                            \"\\n%s%s[%d]:\\n%s%s%s\"
                            % (
                                indent,
                                (\"  \" * (_depth)),
                                i,
                                indent,
                                (\"  \" * (_depth + 1)),
                                str(vv),
                            )
                        )

        return \"\".join(out)

    def pprint(self, *args, **kwargs):
        \"\"\"
        Pretty-printer for parsed results as a list, using the
        `pprint <https://docs.python.org/3/library/pprint.html>`_ module.
        Accepts additional positional or keyword args as defined for
        `pprint.pprint <https://docs.python.org/3/library/pprint.html#pprint.pprint>`_ .

        Example::

            ident = Word(alphas, alphanums)
            num = Word(nums)
            func = Forward()
            term = ident | num | Group('(' + func + ')')
            func <<= ident + Group(Optional(delimited_list(term)))
            result = func.parse_string(\"fna a,b,(fnb c,d,200),100\")
            result.pprint(width=40)

        prints::

            ['fna',
             ['a',
              'b',
              ['(', 'fnb', ['c', 'd', '200'], ')'],
              '100']]
        \"\"\"
        pprint.pprint(self.as_list(), *args, **kwargs)

    # add support for pickle protocol
    def __getstate__(self):
        return (
            self._toklist,
            (
                self._tokdict.copy(),
                self._parent is not None and self._parent() or None,
                self._all_names,
                self._name,
            ),
        )

    def __setstate__(self, state):
        self._toklist, (self._tokdict, par, inAccumNames, self._name) = state
        self._all_names = set(inAccumNames)
        if par is not None:
            self._parent = wkref(par)
        else:
            self._parent = None

    def __getnewargs__(self):
        return self._toklist, self._name

    def __dir__(self):
        return dir(type(self)) + list(self.keys())

    @classmethod
    def from_dict(cls, other, name=None) -> \"ParseResults\":
        \"\"\"
        Helper classmethod to construct a ``ParseResults`` from a ``dict``, preserving the
        name-value relations as results names. If an optional ``name`` argument is
        given, a nested ``ParseResults`` will be returned.
        \"\"\"

        def is_iterable(obj):
            try:
                iter(obj)
            except Exception:
                return False
            else:
                return not isinstance(obj, str_type)

        ret = cls([])
        for k, v in other.items():
            if isinstance(v, Mapping):
                ret += cls.from_dict(v, name=k)
            else:
                ret += cls([v], name=k, asList=is_iterable(v))
        if name is not None:
            ret = cls([ret], name=name)
        return ret

    asList = as_list
    asDict = as_dict
    getName = get_name


MutableMapping.register(ParseResults)
MutableSequence.register(ParseResults)

"""
module_dict["pyparsing/__init__.py"]="""
# module pyparsing.py
#
# Copyright (c) 2003-2022  Paul T. McGuire
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# \"Software\"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__doc__ = \"\"\"
pyparsing module - Classes and methods to define and execute parsing grammars
=============================================================================

The pyparsing module is an alternative approach to creating and
executing simple grammars, vs. the traditional lex/yacc approach, or the
use of regular expressions.  With pyparsing, you don't need to learn
a new syntax for defining grammars or matching expressions - the parsing
module provides a library of classes that you use to construct the
grammar directly in Python.

Here is a program to parse \"Hello, World!\" (or any greeting of the form
``\"<salutation>, <addressee>!\"``), built up using :class:`Word`,
:class:`Literal`, and :class:`And` elements
(the :meth:`'+'<ParserElement.__add__>` operators create :class:`And` expressions,
and the strings are auto-converted to :class:`Literal` expressions)::

    from pyparsing import Word, alphas

    # define grammar of a greeting
    greet = Word(alphas) + \",\" + Word(alphas) + \"!\"

    hello = \"Hello, World!\"
    print(hello, \"->\", greet.parse_string(hello))

The program outputs the following::

    Hello, World! -> ['Hello', ',', 'World', '!']

The Python representation of the grammar is quite readable, owing to the
self-explanatory class names, and the use of :class:`'+'<And>`,
:class:`'|'<MatchFirst>`, :class:`'^'<Or>` and :class:`'&'<Each>` operators.

The :class:`ParseResults` object returned from
:class:`ParserElement.parseString` can be
accessed as a nested list, a dictionary, or an object with named
attributes.

The pyparsing module handles some of the problems that are typically
vexing when writing text parsers:

  - extra or missing whitespace (the above program will also handle
    \"Hello,World!\", \"Hello  ,  World  !\", etc.)
  - quoted strings
  - embedded comments


Getting Started -
-----------------
Visit the classes :class:`ParserElement` and :class:`ParseResults` to
see the base classes that most other pyparsing
classes inherit from. Use the docstrings for examples of how to:

 - construct literal match expressions from :class:`Literal` and
   :class:`CaselessLiteral` classes
 - construct character word-group expressions using the :class:`Word`
   class
 - see how to create repetitive expressions using :class:`ZeroOrMore`
   and :class:`OneOrMore` classes
 - use :class:`'+'<And>`, :class:`'|'<MatchFirst>`, :class:`'^'<Or>`,
   and :class:`'&'<Each>` operators to combine simple expressions into
   more complex ones
 - associate names with your parsed results using
   :class:`ParserElement.setResultsName`
 - access the parsed data, which is returned as a :class:`ParseResults`
   object
 - find some helpful expression short-cuts like :class:`delimitedList`
   and :class:`oneOf`
 - find more useful common expressions in the :class:`pyparsing_common`
   namespace class
\"\"\"
from typing import NamedTuple


class version_info(NamedTuple):
    major: int
    minor: int
    micro: int
    releaselevel: str
    serial: int

    @property
    def __version__(self):
        return (
            \"{}.{}.{}\".format(self.major, self.minor, self.micro)
            + (
                \"{}{}{}\".format(
                    \"r\" if self.releaselevel[0] == \"c\" else \"\",
                    self.releaselevel[0],
                    self.serial,
                ),
                \"\",
            )[self.releaselevel == \"final\"]
        )

    def __str__(self):
        return \"{} {} / {}\".format(__name__, self.__version__, __version_time__)

    def __repr__(self):
        return \"{}.{}({})\".format(
            __name__,
            type(self).__name__,
            \", \".join(\"{}={!r}\".format(*nv) for nv in zip(self._fields, self)),
        )


__version_info__ = version_info(3, 0, 9, \"final\", 0)
__version_time__ = \"05 May 2022 07:02 UTC\"
__version__ = __version_info__.__version__
__versionTime__ = __version_time__
__author__ = \"Paul McGuire <ptmcg.gm+pyparsing@gmail.com>\"

from .util import *
from .exceptions import *
from .actions import *
from .core import __diag__, __compat__
from .results import *
from .core import *
from .core import _builtin_exprs as core_builtin_exprs
from .helpers import *
from .helpers import _builtin_exprs as helper_builtin_exprs

from .unicode import unicode_set, UnicodeRangeList, pyparsing_unicode as unicode
from .testing import pyparsing_test as testing
from .common import (
    pyparsing_common as common,
    _builtin_exprs as common_builtin_exprs,
)

# define backward compat synonyms
if \"pyparsing_unicode\" not in globals():
    pyparsing_unicode = unicode
if \"pyparsing_common\" not in globals():
    pyparsing_common = common
if \"pyparsing_test\" not in globals():
    pyparsing_test = testing

core_builtin_exprs += common_builtin_exprs + helper_builtin_exprs


__all__ = [
    \"__version__\",
    \"__version_time__\",
    \"__author__\",
    \"__compat__\",
    \"__diag__\",
    \"And\",
    \"AtLineStart\",
    \"AtStringStart\",
    \"CaselessKeyword\",
    \"CaselessLiteral\",
    \"CharsNotIn\",
    \"Combine\",
    \"Dict\",
    \"Each\",
    \"Empty\",
    \"FollowedBy\",
    \"Forward\",
    \"GoToColumn\",
    \"Group\",
    \"IndentedBlock\",
    \"Keyword\",
    \"LineEnd\",
    \"LineStart\",
    \"Literal\",
    \"Located\",
    \"PrecededBy\",
    \"MatchFirst\",
    \"NoMatch\",
    \"NotAny\",
    \"OneOrMore\",
    \"OnlyOnce\",
    \"OpAssoc\",
    \"Opt\",
    \"Optional\",
    \"Or\",
    \"ParseBaseException\",
    \"ParseElementEnhance\",
    \"ParseException\",
    \"ParseExpression\",
    \"ParseFatalException\",
    \"ParseResults\",
    \"ParseSyntaxException\",
    \"ParserElement\",
    \"PositionToken\",
    \"QuotedString\",
    \"RecursiveGrammarException\",
    \"Regex\",
    \"SkipTo\",
    \"StringEnd\",
    \"StringStart\",
    \"Suppress\",
    \"Token\",
    \"TokenConverter\",
    \"White\",
    \"Word\",
    \"WordEnd\",
    \"WordStart\",
    \"ZeroOrMore\",
    \"Char\",
    \"alphanums\",
    \"alphas\",
    \"alphas8bit\",
    \"any_close_tag\",
    \"any_open_tag\",
    \"c_style_comment\",
    \"col\",
    \"common_html_entity\",
    \"counted_array\",
    \"cpp_style_comment\",
    \"dbl_quoted_string\",
    \"dbl_slash_comment\",
    \"delimited_list\",
    \"dict_of\",
    \"empty\",
    \"hexnums\",
    \"html_comment\",
    \"identchars\",
    \"identbodychars\",
    \"java_style_comment\",
    \"line\",
    \"line_end\",
    \"line_start\",
    \"lineno\",
    \"make_html_tags\",
    \"make_xml_tags\",
    \"match_only_at_col\",
    \"match_previous_expr\",
    \"match_previous_literal\",
    \"nested_expr\",
    \"null_debug_action\",
    \"nums\",
    \"one_of\",
    \"printables\",
    \"punc8bit\",
    \"python_style_comment\",
    \"quoted_string\",
    \"remove_quotes\",
    \"replace_with\",
    \"replace_html_entity\",
    \"rest_of_line\",
    \"sgl_quoted_string\",
    \"srange\",
    \"string_end\",
    \"string_start\",
    \"trace_parse_action\",
    \"unicode_string\",
    \"with_attribute\",
    \"indentedBlock\",
    \"original_text_for\",
    \"ungroup\",
    \"infix_notation\",
    \"locatedExpr\",
    \"with_class\",
    \"CloseMatch\",
    \"token_map\",
    \"pyparsing_common\",
    \"pyparsing_unicode\",
    \"unicode_set\",
    \"condition_as_parse_action\",
    \"pyparsing_test\",
    # pre-PEP8 compatibility names
    \"__versionTime__\",
    \"anyCloseTag\",
    \"anyOpenTag\",
    \"cStyleComment\",
    \"commonHTMLEntity\",
    \"countedArray\",
    \"cppStyleComment\",
    \"dblQuotedString\",
    \"dblSlashComment\",
    \"delimitedList\",
    \"dictOf\",
    \"htmlComment\",
    \"javaStyleComment\",
    \"lineEnd\",
    \"lineStart\",
    \"makeHTMLTags\",
    \"makeXMLTags\",
    \"matchOnlyAtCol\",
    \"matchPreviousExpr\",
    \"matchPreviousLiteral\",
    \"nestedExpr\",
    \"nullDebugAction\",
    \"oneOf\",
    \"opAssoc\",
    \"pythonStyleComment\",
    \"quotedString\",
    \"removeQuotes\",
    \"replaceHTMLEntity\",
    \"replaceWith\",
    \"restOfLine\",
    \"sglQuotedString\",
    \"stringEnd\",
    \"stringStart\",
    \"traceParseAction\",
    \"unicodeString\",
    \"withAttribute\",
    \"indentedBlock\",
    \"originalTextFor\",
    \"infixNotation\",
    \"locatedExpr\",
    \"withClass\",
    \"tokenMap\",
    \"conditionAsParseAction\",
    \"autoname_elements\",
]

"""
module_dict["pyparsing/util.py"]="""
# util.py
import warnings
import types
import collections
import itertools
from functools import lru_cache
from typing import List, Union, Iterable

_bslash = chr(92)


class __config_flags:
    \"\"\"Internal class for defining compatibility and debugging flags\"\"\"

    _all_names: List[str] = []
    _fixed_names: List[str] = []
    _type_desc = \"configuration\"

    @classmethod
    def _set(cls, dname, value):
        if dname in cls._fixed_names:
            warnings.warn(
                \"{}.{} {} is {} and cannot be overridden\".format(
                    cls.__name__,
                    dname,
                    cls._type_desc,
                    str(getattr(cls, dname)).upper(),
                )
            )
            return
        if dname in cls._all_names:
            setattr(cls, dname, value)
        else:
            raise ValueError(\"no such {} {!r}\".format(cls._type_desc, dname))

    enable = classmethod(lambda cls, name: cls._set(name, True))
    disable = classmethod(lambda cls, name: cls._set(name, False))


@lru_cache(maxsize=128)
def col(loc: int, strg: str) -> int:
    \"\"\"
    Returns current column within a string, counting newlines as line separators.
    The first column is number 1.

    Note: the default parsing behavior is to expand tabs in the input string
    before starting the parsing process.  See
    :class:`ParserElement.parseString` for more
    information on parsing strings containing ``<TAB>`` s, and suggested
    methods to maintain a consistent view of the parsed string, the parse
    location, and line and column positions within the parsed string.
    \"\"\"
    s = strg
    return 1 if 0 < loc < len(s) and s[loc - 1] == \"\\n\" else loc - s.rfind(\"\\n\", 0, loc)


@lru_cache(maxsize=128)
def lineno(loc: int, strg: str) -> int:
    \"\"\"Returns current line number within a string, counting newlines as line separators.
    The first line is number 1.

    Note - the default parsing behavior is to expand tabs in the input string
    before starting the parsing process.  See :class:`ParserElement.parseString`
    for more information on parsing strings containing ``<TAB>`` s, and
    suggested methods to maintain a consistent view of the parsed string, the
    parse location, and line and column positions within the parsed string.
    \"\"\"
    return strg.count(\"\\n\", 0, loc) + 1


@lru_cache(maxsize=128)
def line(loc: int, strg: str) -> str:
    \"\"\"
    Returns the line of text containing loc within a string, counting newlines as line separators.
    \"\"\"
    last_cr = strg.rfind(\"\\n\", 0, loc)
    next_cr = strg.find(\"\\n\", loc)
    return strg[last_cr + 1 : next_cr] if next_cr >= 0 else strg[last_cr + 1 :]


class _UnboundedCache:
    def __init__(self):
        cache = {}
        cache_get = cache.get
        self.not_in_cache = not_in_cache = object()

        def get(_, key):
            return cache_get(key, not_in_cache)

        def set_(_, key, value):
            cache[key] = value

        def clear(_):
            cache.clear()

        self.size = None
        self.get = types.MethodType(get, self)
        self.set = types.MethodType(set_, self)
        self.clear = types.MethodType(clear, self)


class _FifoCache:
    def __init__(self, size):
        self.not_in_cache = not_in_cache = object()
        cache = collections.OrderedDict()
        cache_get = cache.get

        def get(_, key):
            return cache_get(key, not_in_cache)

        def set_(_, key, value):
            cache[key] = value
            while len(cache) > size:
                cache.popitem(last=False)

        def clear(_):
            cache.clear()

        self.size = size
        self.get = types.MethodType(get, self)
        self.set = types.MethodType(set_, self)
        self.clear = types.MethodType(clear, self)


class LRUMemo:
    \"\"\"
    A memoizing mapping that retains `capacity` deleted items

    The memo tracks retained items by their access order; once `capacity` items
    are retained, the least recently used item is discarded.
    \"\"\"

    def __init__(self, capacity):
        self._capacity = capacity
        self._active = {}
        self._memory = collections.OrderedDict()

    def __getitem__(self, key):
        try:
            return self._active[key]
        except KeyError:
            self._memory.move_to_end(key)
            return self._memory[key]

    def __setitem__(self, key, value):
        self._memory.pop(key, None)
        self._active[key] = value

    def __delitem__(self, key):
        try:
            value = self._active.pop(key)
        except KeyError:
            pass
        else:
            while len(self._memory) >= self._capacity:
                self._memory.popitem(last=False)
            self._memory[key] = value

    def clear(self):
        self._active.clear()
        self._memory.clear()


class UnboundedMemo(dict):
    \"\"\"
    A memoizing mapping that retains all deleted items
    \"\"\"

    def __delitem__(self, key):
        pass


def _escape_regex_range_chars(s: str) -> str:
    # escape these chars: ^-[]
    for c in r\"\\^-[]\":
        s = s.replace(c, _bslash + c)
    s = s.replace(\"\\n\", r\"\\n\")
    s = s.replace(\"\\t\", r\"\\t\")
    return str(s)


def _collapse_string_to_ranges(
    s: Union[str, Iterable[str]], re_escape: bool = True
) -> str:
    def is_consecutive(c):
        c_int = ord(c)
        is_consecutive.prev, prev = c_int, is_consecutive.prev
        if c_int - prev > 1:
            is_consecutive.value = next(is_consecutive.counter)
        return is_consecutive.value

    is_consecutive.prev = 0
    is_consecutive.counter = itertools.count()
    is_consecutive.value = -1

    def escape_re_range_char(c):
        return \"\\\\\" + c if c in r\"\\^-][\" else c

    def no_escape_re_range_char(c):
        return c

    if not re_escape:
        escape_re_range_char = no_escape_re_range_char

    ret = []
    s = \"\".join(sorted(set(s)))
    if len(s) > 3:
        for _, chars in itertools.groupby(s, key=is_consecutive):
            first = last = next(chars)
            last = collections.deque(
                itertools.chain(iter([last]), chars), maxlen=1
            ).pop()
            if first == last:
                ret.append(escape_re_range_char(first))
            else:
                sep = \"\" if ord(last) == ord(first) + 1 else \"-\"
                ret.append(
                    \"{}{}{}\".format(
                        escape_re_range_char(first), sep, escape_re_range_char(last)
                    )
                )
    else:
        ret = [escape_re_range_char(c) for c in s]

    return \"\".join(ret)


def _flatten(ll: list) -> list:
    ret = []
    for i in ll:
        if isinstance(i, list):
            ret.extend(_flatten(i))
        else:
            ret.append(i)
    return ret

"""
module_dict["pyparsing/common.py"]="""
# common.py
from .core import *
from .helpers import delimited_list, any_open_tag, any_close_tag
from datetime import datetime


# some other useful expressions - using lower-case class name since we are really using this as a namespace
class pyparsing_common:
    \"\"\"Here are some common low-level expressions that may be useful in
    jump-starting parser development:

    - numeric forms (:class:`integers<integer>`, :class:`reals<real>`,
      :class:`scientific notation<sci_real>`)
    - common :class:`programming identifiers<identifier>`
    - network addresses (:class:`MAC<mac_address>`,
      :class:`IPv4<ipv4_address>`, :class:`IPv6<ipv6_address>`)
    - ISO8601 :class:`dates<iso8601_date>` and
      :class:`datetime<iso8601_datetime>`
    - :class:`UUID<uuid>`
    - :class:`comma-separated list<comma_separated_list>`
    - :class:`url`

    Parse actions:

    - :class:`convertToInteger`
    - :class:`convertToFloat`
    - :class:`convertToDate`
    - :class:`convertToDatetime`
    - :class:`stripHTMLTags`
    - :class:`upcaseTokens`
    - :class:`downcaseTokens`

    Example::

        pyparsing_common.number.runTests('''
            # any int or real number, returned as the appropriate type
            100
            -100
            +100
            3.14159
            6.02e23
            1e-12
            ''')

        pyparsing_common.fnumber.runTests('''
            # any int or real number, returned as float
            100
            -100
            +100
            3.14159
            6.02e23
            1e-12
            ''')

        pyparsing_common.hex_integer.runTests('''
            # hex numbers
            100
            FF
            ''')

        pyparsing_common.fraction.runTests('''
            # fractions
            1/2
            -3/4
            ''')

        pyparsing_common.mixed_integer.runTests('''
            # mixed fractions
            1
            1/2
            -3/4
            1-3/4
            ''')

        import uuid
        pyparsing_common.uuid.setParseAction(tokenMap(uuid.UUID))
        pyparsing_common.uuid.runTests('''
            # uuid
            12345678-1234-5678-1234-567812345678
            ''')

    prints::

        # any int or real number, returned as the appropriate type
        100
        [100]

        -100
        [-100]

        +100
        [100]

        3.14159
        [3.14159]

        6.02e23
        [6.02e+23]

        1e-12
        [1e-12]

        # any int or real number, returned as float
        100
        [100.0]

        -100
        [-100.0]

        +100
        [100.0]

        3.14159
        [3.14159]

        6.02e23
        [6.02e+23]

        1e-12
        [1e-12]

        # hex numbers
        100
        [256]

        FF
        [255]

        # fractions
        1/2
        [0.5]

        -3/4
        [-0.75]

        # mixed fractions
        1
        [1]

        1/2
        [0.5]

        -3/4
        [-0.75]

        1-3/4
        [1.75]

        # uuid
        12345678-1234-5678-1234-567812345678
        [UUID('12345678-1234-5678-1234-567812345678')]
    \"\"\"

    convert_to_integer = token_map(int)
    \"\"\"
    Parse action for converting parsed integers to Python int
    \"\"\"

    convert_to_float = token_map(float)
    \"\"\"
    Parse action for converting parsed numbers to Python float
    \"\"\"

    integer = Word(nums).set_name(\"integer\").set_parse_action(convert_to_integer)
    \"\"\"expression that parses an unsigned integer, returns an int\"\"\"

    hex_integer = (
        Word(hexnums).set_name(\"hex integer\").set_parse_action(token_map(int, 16))
    )
    \"\"\"expression that parses a hexadecimal integer, returns an int\"\"\"

    signed_integer = (
        Regex(r\"[+-]?\\d+\")
        .set_name(\"signed integer\")
        .set_parse_action(convert_to_integer)
    )
    \"\"\"expression that parses an integer with optional leading sign, returns an int\"\"\"

    fraction = (
        signed_integer().set_parse_action(convert_to_float)
        + \"/\"
        + signed_integer().set_parse_action(convert_to_float)
    ).set_name(\"fraction\")
    \"\"\"fractional expression of an integer divided by an integer, returns a float\"\"\"
    fraction.add_parse_action(lambda tt: tt[0] / tt[-1])

    mixed_integer = (
        fraction | signed_integer + Opt(Opt(\"-\").suppress() + fraction)
    ).set_name(\"fraction or mixed integer-fraction\")
    \"\"\"mixed integer of the form 'integer - fraction', with optional leading integer, returns float\"\"\"
    mixed_integer.add_parse_action(sum)

    real = (
        Regex(r\"[+-]?(?:\\d+\\.\\d*|\\.\\d+)\")
        .set_name(\"real number\")
        .set_parse_action(convert_to_float)
    )
    \"\"\"expression that parses a floating point number and returns a float\"\"\"

    sci_real = (
        Regex(r\"[+-]?(?:\\d+(?:[eE][+-]?\\d+)|(?:\\d+\\.\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?)\")
        .set_name(\"real number with scientific notation\")
        .set_parse_action(convert_to_float)
    )
    \"\"\"expression that parses a floating point number with optional
    scientific notation and returns a float\"\"\"

    # streamlining this expression makes the docs nicer-looking
    number = (sci_real | real | signed_integer).setName(\"number\").streamline()
    \"\"\"any numeric expression, returns the corresponding Python type\"\"\"

    fnumber = (
        Regex(r\"[+-]?\\d+\\.?\\d*([eE][+-]?\\d+)?\")
        .set_name(\"fnumber\")
        .set_parse_action(convert_to_float)
    )
    \"\"\"any int or real number, returned as float\"\"\"

    identifier = Word(identchars, identbodychars).set_name(\"identifier\")
    \"\"\"typical code identifier (leading alpha or '_', followed by 0 or more alphas, nums, or '_')\"\"\"

    ipv4_address = Regex(
        r\"(25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})(\\.(25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})){3}\"
    ).set_name(\"IPv4 address\")
    \"IPv4 address (``0.0.0.0 - 255.255.255.255``)\"

    _ipv6_part = Regex(r\"[0-9a-fA-F]{1,4}\").set_name(\"hex_integer\")
    _full_ipv6_address = (_ipv6_part + (\":\" + _ipv6_part) * 7).set_name(
        \"full IPv6 address\"
    )
    _short_ipv6_address = (
        Opt(_ipv6_part + (\":\" + _ipv6_part) * (0, 6))
        + \"::\"
        + Opt(_ipv6_part + (\":\" + _ipv6_part) * (0, 6))
    ).set_name(\"short IPv6 address\")
    _short_ipv6_address.add_condition(
        lambda t: sum(1 for tt in t if pyparsing_common._ipv6_part.matches(tt)) < 8
    )
    _mixed_ipv6_address = (\"::ffff:\" + ipv4_address).set_name(\"mixed IPv6 address\")
    ipv6_address = Combine(
        (_full_ipv6_address | _mixed_ipv6_address | _short_ipv6_address).set_name(
            \"IPv6 address\"
        )
    ).set_name(\"IPv6 address\")
    \"IPv6 address (long, short, or mixed form)\"

    mac_address = Regex(
        r\"[0-9a-fA-F]{2}([:.-])[0-9a-fA-F]{2}(?:\\1[0-9a-fA-F]{2}){4}\"
    ).set_name(\"MAC address\")
    \"MAC address xx:xx:xx:xx:xx (may also have '-' or '.' delimiters)\"

    @staticmethod
    def convert_to_date(fmt: str = \"%Y-%m-%d\"):
        \"\"\"
        Helper to create a parse action for converting parsed date string to Python datetime.date

        Params -
        - fmt - format to be passed to datetime.strptime (default= ``\"%Y-%m-%d\"``)

        Example::

            date_expr = pyparsing_common.iso8601_date.copy()
            date_expr.setParseAction(pyparsing_common.convertToDate())
            print(date_expr.parseString(\"1999-12-31\"))

        prints::

            [datetime.date(1999, 12, 31)]
        \"\"\"

        def cvt_fn(ss, ll, tt):
            try:
                return datetime.strptime(tt[0], fmt).date()
            except ValueError as ve:
                raise ParseException(ss, ll, str(ve))

        return cvt_fn

    @staticmethod
    def convert_to_datetime(fmt: str = \"%Y-%m-%dT%H:%M:%S.%f\"):
        \"\"\"Helper to create a parse action for converting parsed
        datetime string to Python datetime.datetime

        Params -
        - fmt - format to be passed to datetime.strptime (default= ``\"%Y-%m-%dT%H:%M:%S.%f\"``)

        Example::

            dt_expr = pyparsing_common.iso8601_datetime.copy()
            dt_expr.setParseAction(pyparsing_common.convertToDatetime())
            print(dt_expr.parseString(\"1999-12-31T23:59:59.999\"))

        prints::

            [datetime.datetime(1999, 12, 31, 23, 59, 59, 999000)]
        \"\"\"

        def cvt_fn(s, l, t):
            try:
                return datetime.strptime(t[0], fmt)
            except ValueError as ve:
                raise ParseException(s, l, str(ve))

        return cvt_fn

    iso8601_date = Regex(
        r\"(?P<year>\\d{4})(?:-(?P<month>\\d\\d)(?:-(?P<day>\\d\\d))?)?\"
    ).set_name(\"ISO8601 date\")
    \"ISO8601 date (``yyyy-mm-dd``)\"

    iso8601_datetime = Regex(
        r\"(?P<year>\\d{4})-(?P<month>\\d\\d)-(?P<day>\\d\\d)[T ](?P<hour>\\d\\d):(?P<minute>\\d\\d)(:(?P<second>\\d\\d(\\.\\d*)?)?)?(?P<tz>Z|[+-]\\d\\d:?\\d\\d)?\"
    ).set_name(\"ISO8601 datetime\")
    \"ISO8601 datetime (``yyyy-mm-ddThh:mm:ss.s(Z|+-00:00)``) - trailing seconds, milliseconds, and timezone optional; accepts separating ``'T'`` or ``' '``\"

    uuid = Regex(r\"[0-9a-fA-F]{8}(-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}\").set_name(\"UUID\")
    \"UUID (``xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx``)\"

    _html_stripper = any_open_tag.suppress() | any_close_tag.suppress()

    @staticmethod
    def strip_html_tags(s: str, l: int, tokens: ParseResults):
        \"\"\"Parse action to remove HTML tags from web page HTML source

        Example::

            # strip HTML links from normal text
            text = '<td>More info at the <a href=\"https://github.com/pyparsing/pyparsing/wiki\">pyparsing</a> wiki page</td>'
            td, td_end = makeHTMLTags(\"TD\")
            table_text = td + SkipTo(td_end).setParseAction(pyparsing_common.stripHTMLTags)(\"body\") + td_end
            print(table_text.parseString(text).body)

        Prints::

            More info at the pyparsing wiki page
        \"\"\"
        return pyparsing_common._html_stripper.transform_string(tokens[0])

    _commasepitem = (
        Combine(
            OneOrMore(
                ~Literal(\",\")
                + ~LineEnd()
                + Word(printables, exclude_chars=\",\")
                + Opt(White(\" \\t\") + ~FollowedBy(LineEnd() | \",\"))
            )
        )
        .streamline()
        .set_name(\"commaItem\")
    )
    comma_separated_list = delimited_list(
        Opt(quoted_string.copy() | _commasepitem, default=\"\")
    ).set_name(\"comma separated list\")
    \"\"\"Predefined expression of 1 or more printable words or quoted strings, separated by commas.\"\"\"

    upcase_tokens = staticmethod(token_map(lambda t: t.upper()))
    \"\"\"Parse action to convert tokens to upper case.\"\"\"

    downcase_tokens = staticmethod(token_map(lambda t: t.lower()))
    \"\"\"Parse action to convert tokens to lower case.\"\"\"

    # fmt: off
    url = Regex(
        # https://mathiasbynens.be/demo/url-regex
        # https://gist.github.com/dperini/729294
        r\"^\" +
        # protocol identifier (optional)
        # short syntax // still required
        r\"(?:(?:(?P<scheme>https?|ftp):)?\\/\\/)\" +
        # user:pass BasicAuth (optional)
        r\"(?:(?P<auth>\\S+(?::\\S*)?)@)?\" +
        r\"(?P<host>\" +
        # IP address exclusion
        # private & local networks
        r\"(?!(?:10|127)(?:\\.\\d{1,3}){3})\" +
        r\"(?!(?:169\\.254|192\\.168)(?:\\.\\d{1,3}){2})\" +
        r\"(?!172\\.(?:1[6-9]|2\\d|3[0-1])(?:\\.\\d{1,3}){2})\" +
        # IP address dotted notation octets
        # excludes loopback network 0.0.0.0
        # excludes reserved space >= 224.0.0.0
        # excludes network & broadcast addresses
        # (first & last IP address of each class)
        r\"(?:[1-9]\\d?|1\\d\\d|2[01]\\d|22[0-3])\" +
        r\"(?:\\.(?:1?\\d{1,2}|2[0-4]\\d|25[0-5])){2}\" +
        r\"(?:\\.(?:[1-9]\\d?|1\\d\\d|2[0-4]\\d|25[0-4]))\" +
        r\"|\" +
        # host & domain names, may end with dot
        # can be replaced by a shortest alternative
        # (?![-_])(?:[-\\w\\u00a1-\\uffff]{0,63}[^-_]\\.)+
        r\"(?:\" +
        r\"(?:\" +
        r\"[a-z0-9\\u00a1-\\uffff]\" +
        r\"[a-z0-9\\u00a1-\\uffff_-]{0,62}\" +
        r\")?\" +
        r\"[a-z0-9\\u00a1-\\uffff]\\.\" +
        r\")+\" +
        # TLD identifier name, may end with dot
        r\"(?:[a-z\\u00a1-\\uffff]{2,}\\.?)\" +
        r\")\" +
        # port number (optional)
        r\"(:(?P<port>\\d{2,5}))?\" +
        # resource path (optional)
        r\"(?P<path>\\/[^?# ]*)?\" +
        # query string (optional)
        r\"(\\?(?P<query>[^#]*))?\" +
        # fragment (optional)
        r\"(#(?P<fragment>\\S*))?\" +
        r\"$\"
    ).set_name(\"url\")
    # fmt: on

    # pre-PEP8 compatibility names
    convertToInteger = convert_to_integer
    convertToFloat = convert_to_float
    convertToDate = convert_to_date
    convertToDatetime = convert_to_datetime
    stripHTMLTags = strip_html_tags
    upcaseTokens = upcase_tokens
    downcaseTokens = downcase_tokens


_builtin_exprs = [
    v for v in vars(pyparsing_common).values() if isinstance(v, ParserElement)
]

"""
module_dict["pyparsing/actions.py"]="""
# actions.py

from .exceptions import ParseException
from .util import col


class OnlyOnce:
    \"\"\"
    Wrapper for parse actions, to ensure they are only called once.
    \"\"\"

    def __init__(self, method_call):
        from .core import _trim_arity

        self.callable = _trim_arity(method_call)
        self.called = False

    def __call__(self, s, l, t):
        if not self.called:
            results = self.callable(s, l, t)
            self.called = True
            return results
        raise ParseException(s, l, \"OnlyOnce obj called multiple times w/out reset\")

    def reset(self):
        \"\"\"
        Allow the associated parse action to be called once more.
        \"\"\"

        self.called = False


def match_only_at_col(n):
    \"\"\"
    Helper method for defining parse actions that require matching at
    a specific column in the input text.
    \"\"\"

    def verify_col(strg, locn, toks):
        if col(locn, strg) != n:
            raise ParseException(strg, locn, \"matched token not at column {}\".format(n))

    return verify_col


def replace_with(repl_str):
    \"\"\"
    Helper method for common parse actions that simply return
    a literal value.  Especially useful when used with
    :class:`transform_string<ParserElement.transform_string>` ().

    Example::

        num = Word(nums).set_parse_action(lambda toks: int(toks[0]))
        na = one_of(\"N/A NA\").set_parse_action(replace_with(math.nan))
        term = na | num

        term[1, ...].parse_string(\"324 234 N/A 234\") # -> [324, 234, nan, 234]
    \"\"\"
    return lambda s, l, t: [repl_str]


def remove_quotes(s, l, t):
    \"\"\"
    Helper parse action for removing quotation marks from parsed
    quoted strings.

    Example::

        # by default, quotation marks are included in parsed results
        quoted_string.parse_string(\"'Now is the Winter of our Discontent'\") # -> [\"'Now is the Winter of our Discontent'\"]

        # use remove_quotes to strip quotation marks from parsed results
        quoted_string.set_parse_action(remove_quotes)
        quoted_string.parse_string(\"'Now is the Winter of our Discontent'\") # -> [\"Now is the Winter of our Discontent\"]
    \"\"\"
    return t[0][1:-1]


def with_attribute(*args, **attr_dict):
    \"\"\"
    Helper to create a validating parse action to be used with start
    tags created with :class:`make_xml_tags` or
    :class:`make_html_tags`. Use ``with_attribute`` to qualify
    a starting tag with a required attribute value, to avoid false
    matches on common tags such as ``<TD>`` or ``<DIV>``.

    Call ``with_attribute`` with a series of attribute names and
    values. Specify the list of filter attributes names and values as:

    - keyword arguments, as in ``(align=\"right\")``, or
    - as an explicit dict with ``**`` operator, when an attribute
      name is also a Python reserved word, as in ``**{\"class\":\"Customer\", \"align\":\"right\"}``
    - a list of name-value tuples, as in ``((\"ns1:class\", \"Customer\"), (\"ns2:align\", \"right\"))``

    For attribute names with a namespace prefix, you must use the second
    form.  Attribute names are matched insensitive to upper/lower case.

    If just testing for ``class`` (with or without a namespace), use
    :class:`with_class`.

    To verify that the attribute exists, but without specifying a value,
    pass ``with_attribute.ANY_VALUE`` as the value.

    Example::

        html = '''
            <div>
            Some text
            <div type=\"grid\">1 4 0 1 0</div>
            <div type=\"graph\">1,3 2,3 1,1</div>
            <div>this has no type</div>
            </div>

        '''
        div,div_end = make_html_tags(\"div\")

        # only match div tag having a type attribute with value \"grid\"
        div_grid = div().set_parse_action(with_attribute(type=\"grid\"))
        grid_expr = div_grid + SkipTo(div | div_end)(\"body\")
        for grid_header in grid_expr.search_string(html):
            print(grid_header.body)

        # construct a match with any div tag having a type attribute, regardless of the value
        div_any_type = div().set_parse_action(with_attribute(type=with_attribute.ANY_VALUE))
        div_expr = div_any_type + SkipTo(div | div_end)(\"body\")
        for div_header in div_expr.search_string(html):
            print(div_header.body)

    prints::

        1 4 0 1 0

        1 4 0 1 0
        1,3 2,3 1,1
    \"\"\"
    if args:
        attrs = args[:]
    else:
        attrs = attr_dict.items()
    attrs = [(k, v) for k, v in attrs]

    def pa(s, l, tokens):
        for attrName, attrValue in attrs:
            if attrName not in tokens:
                raise ParseException(s, l, \"no matching attribute \" + attrName)
            if attrValue != with_attribute.ANY_VALUE and tokens[attrName] != attrValue:
                raise ParseException(
                    s,
                    l,
                    \"attribute {!r} has value {!r}, must be {!r}\".format(
                        attrName, tokens[attrName], attrValue
                    ),
                )

    return pa


with_attribute.ANY_VALUE = object()


def with_class(classname, namespace=\"\"):
    \"\"\"
    Simplified version of :class:`with_attribute` when
    matching on a div class - made difficult because ``class`` is
    a reserved word in Python.

    Example::

        html = '''
            <div>
            Some text
            <div class=\"grid\">1 4 0 1 0</div>
            <div class=\"graph\">1,3 2,3 1,1</div>
            <div>this &lt;div&gt; has no class</div>
            </div>

        '''
        div,div_end = make_html_tags(\"div\")
        div_grid = div().set_parse_action(with_class(\"grid\"))

        grid_expr = div_grid + SkipTo(div | div_end)(\"body\")
        for grid_header in grid_expr.search_string(html):
            print(grid_header.body)

        div_any_type = div().set_parse_action(with_class(withAttribute.ANY_VALUE))
        div_expr = div_any_type + SkipTo(div | div_end)(\"body\")
        for div_header in div_expr.search_string(html):
            print(div_header.body)

    prints::

        1 4 0 1 0

        1 4 0 1 0
        1,3 2,3 1,1
    \"\"\"
    classattr = \"{}:class\".format(namespace) if namespace else \"class\"
    return with_attribute(**{classattr: classname})


# pre-PEP8 compatibility symbols
replaceWith = replace_with
removeQuotes = remove_quotes
withAttribute = with_attribute
withClass = with_class
matchOnlyAtCol = match_only_at_col

"""
module_dict["pyparsing/diagram/__init__.py"]="""
import railroad
import pyparsing
import typing
from typing import (
    List,
    NamedTuple,
    Generic,
    TypeVar,
    Dict,
    Callable,
    Set,
    Iterable,
)
from jinja2 import Template
from io import StringIO
import inspect


jinja2_template_source = \"\"\"\\
<!DOCTYPE html>
<html>
<head>
    {% if not head %}
        <style type=\"text/css\">
            .railroad-heading {
                font-family: monospace;
            }
        </style>
    {% else %}
        {{ head | safe }}
    {% endif %}
</head>
<body>
{{ body | safe }}
{% for diagram in diagrams %}
    <div class=\"railroad-group\">
        <h1 class=\"railroad-heading\">{{ diagram.title }}</h1>
        <div class=\"railroad-description\">{{ diagram.text }}</div>
        <div class=\"railroad-svg\">
            {{ diagram.svg }}
        </div>
    </div>
{% endfor %}
</body>
</html>
\"\"\"

template = Template(jinja2_template_source)

# Note: ideally this would be a dataclass, but we're supporting Python 3.5+ so we can't do this yet
NamedDiagram = NamedTuple(
    \"NamedDiagram\",
    [(\"name\", str), (\"diagram\", typing.Optional[railroad.DiagramItem]), (\"index\", int)],
)
\"\"\"
A simple structure for associating a name with a railroad diagram
\"\"\"

T = TypeVar(\"T\")


class EachItem(railroad.Group):
    \"\"\"
    Custom railroad item to compose a:
    - Group containing a
      - OneOrMore containing a
        - Choice of the elements in the Each
    with the group label indicating that all must be matched
    \"\"\"

    all_label = \"[ALL]\"

    def __init__(self, *items):
        choice_item = railroad.Choice(len(items) - 1, *items)
        one_or_more_item = railroad.OneOrMore(item=choice_item)
        super().__init__(one_or_more_item, label=self.all_label)


class AnnotatedItem(railroad.Group):
    \"\"\"
    Simple subclass of Group that creates an annotation label
    \"\"\"

    def __init__(self, label: str, item):
        super().__init__(item=item, label=\"[{}]\".format(label) if label else label)


class EditablePartial(Generic[T]):
    \"\"\"
    Acts like a functools.partial, but can be edited. In other words, it represents a type that hasn't yet been
    constructed.
    \"\"\"

    # We need this here because the railroad constructors actually transform the data, so can't be called until the
    # entire tree is assembled

    def __init__(self, func: Callable[..., T], args: list, kwargs: dict):
        self.func = func
        self.args = args
        self.kwargs = kwargs

    @classmethod
    def from_call(cls, func: Callable[..., T], *args, **kwargs) -> \"EditablePartial[T]\":
        \"\"\"
        If you call this function in the same way that you would call the constructor, it will store the arguments
        as you expect. For example EditablePartial.from_call(Fraction, 1, 3)() == Fraction(1, 3)
        \"\"\"
        return EditablePartial(func=func, args=list(args), kwargs=kwargs)

    @property
    def name(self):
        return self.kwargs[\"name\"]

    def __call__(self) -> T:
        \"\"\"
        Evaluate the partial and return the result
        \"\"\"
        args = self.args.copy()
        kwargs = self.kwargs.copy()

        # This is a helpful hack to allow you to specify varargs parameters (e.g. *args) as keyword args (e.g.
        # args=['list', 'of', 'things'])
        arg_spec = inspect.getfullargspec(self.func)
        if arg_spec.varargs in self.kwargs:
            args += kwargs.pop(arg_spec.varargs)

        return self.func(*args, **kwargs)


def railroad_to_html(diagrams: List[NamedDiagram], **kwargs) -> str:
    \"\"\"
    Given a list of NamedDiagram, produce a single HTML string that visualises those diagrams
    :params kwargs: kwargs to be passed in to the template
    \"\"\"
    data = []
    for diagram in diagrams:
        if diagram.diagram is None:
            continue
        io = StringIO()
        diagram.diagram.writeSvg(io.write)
        title = diagram.name
        if diagram.index == 0:
            title += \" (root)\"
        data.append({\"title\": title, \"text\": \"\", \"svg\": io.getvalue()})

    return template.render(diagrams=data, **kwargs)


def resolve_partial(partial: \"EditablePartial[T]\") -> T:
    \"\"\"
    Recursively resolves a collection of Partials into whatever type they are
    \"\"\"
    if isinstance(partial, EditablePartial):
        partial.args = resolve_partial(partial.args)
        partial.kwargs = resolve_partial(partial.kwargs)
        return partial()
    elif isinstance(partial, list):
        return [resolve_partial(x) for x in partial]
    elif isinstance(partial, dict):
        return {key: resolve_partial(x) for key, x in partial.items()}
    else:
        return partial


def to_railroad(
    element: pyparsing.ParserElement,
    diagram_kwargs: typing.Optional[dict] = None,
    vertical: int = 3,
    show_results_names: bool = False,
    show_groups: bool = False,
) -> List[NamedDiagram]:
    \"\"\"
    Convert a pyparsing element tree into a list of diagrams. This is the recommended entrypoint to diagram
    creation if you want to access the Railroad tree before it is converted to HTML
    :param element: base element of the parser being diagrammed
    :param diagram_kwargs: kwargs to pass to the Diagram() constructor
    :param vertical: (optional) - int - limit at which number of alternatives should be
       shown vertically instead of horizontally
    :param show_results_names - bool to indicate whether results name annotations should be
       included in the diagram
    :param show_groups - bool to indicate whether groups should be highlighted with an unlabeled
       surrounding box
    \"\"\"
    # Convert the whole tree underneath the root
    lookup = ConverterState(diagram_kwargs=diagram_kwargs or {})
    _to_diagram_element(
        element,
        lookup=lookup,
        parent=None,
        vertical=vertical,
        show_results_names=show_results_names,
        show_groups=show_groups,
    )

    root_id = id(element)
    # Convert the root if it hasn't been already
    if root_id in lookup:
        if not element.customName:
            lookup[root_id].name = \"\"
        lookup[root_id].mark_for_extraction(root_id, lookup, force=True)

    # Now that we're finished, we can convert from intermediate structures into Railroad elements
    diags = list(lookup.diagrams.values())
    if len(diags) > 1:
        # collapse out duplicate diags with the same name
        seen = set()
        deduped_diags = []
        for d in diags:
            # don't extract SkipTo elements, they are uninformative as subdiagrams
            if d.name == \"...\":
                continue
            if d.name is not None and d.name not in seen:
                seen.add(d.name)
                deduped_diags.append(d)
        resolved = [resolve_partial(partial) for partial in deduped_diags]
    else:
        # special case - if just one diagram, always display it, even if
        # it has no name
        resolved = [resolve_partial(partial) for partial in diags]
    return sorted(resolved, key=lambda diag: diag.index)


def _should_vertical(
    specification: int, exprs: Iterable[pyparsing.ParserElement]
) -> bool:
    \"\"\"
    Returns true if we should return a vertical list of elements
    \"\"\"
    if specification is None:
        return False
    else:
        return len(_visible_exprs(exprs)) >= specification


class ElementState:
    \"\"\"
    State recorded for an individual pyparsing Element
    \"\"\"

    # Note: this should be a dataclass, but we have to support Python 3.5
    def __init__(
        self,
        element: pyparsing.ParserElement,
        converted: EditablePartial,
        parent: EditablePartial,
        number: int,
        name: str = None,
        parent_index: typing.Optional[int] = None,
    ):
        #: The pyparsing element that this represents
        self.element: pyparsing.ParserElement = element
        #: The name of the element
        self.name: typing.Optional[str] = name
        #: The output Railroad element in an unconverted state
        self.converted: EditablePartial = converted
        #: The parent Railroad element, which we store so that we can extract this if it's duplicated
        self.parent: EditablePartial = parent
        #: The order in which we found this element, used for sorting diagrams if this is extracted into a diagram
        self.number: int = number
        #: The index of this inside its parent
        self.parent_index: typing.Optional[int] = parent_index
        #: If true, we should extract this out into a subdiagram
        self.extract: bool = False
        #: If true, all of this element's children have been filled out
        self.complete: bool = False

    def mark_for_extraction(
        self, el_id: int, state: \"ConverterState\", name: str = None, force: bool = False
    ):
        \"\"\"
        Called when this instance has been seen twice, and thus should eventually be extracted into a sub-diagram
        :param el_id: id of the element
        :param state: element/diagram state tracker
        :param name: name to use for this element's text
        :param force: If true, force extraction now, regardless of the state of this. Only useful for extracting the
        root element when we know we're finished
        \"\"\"
        self.extract = True

        # Set the name
        if not self.name:
            if name:
                # Allow forcing a custom name
                self.name = name
            elif self.element.customName:
                self.name = self.element.customName
            else:
                self.name = \"\"

        # Just because this is marked for extraction doesn't mean we can do it yet. We may have to wait for children
        # to be added
        # Also, if this is just a string literal etc, don't bother extracting it
        if force or (self.complete and _worth_extracting(self.element)):
            state.extract_into_diagram(el_id)


class ConverterState:
    \"\"\"
    Stores some state that persists between recursions into the element tree
    \"\"\"

    def __init__(self, diagram_kwargs: typing.Optional[dict] = None):
        #: A dictionary mapping ParserElements to state relating to them
        self._element_diagram_states: Dict[int, ElementState] = {}
        #: A dictionary mapping ParserElement IDs to subdiagrams generated from them
        self.diagrams: Dict[int, EditablePartial[NamedDiagram]] = {}
        #: The index of the next unnamed element
        self.unnamed_index: int = 1
        #: The index of the next element. This is used for sorting
        self.index: int = 0
        #: Shared kwargs that are used to customize the construction of diagrams
        self.diagram_kwargs: dict = diagram_kwargs or {}
        self.extracted_diagram_names: Set[str] = set()

    def __setitem__(self, key: int, value: ElementState):
        self._element_diagram_states[key] = value

    def __getitem__(self, key: int) -> ElementState:
        return self._element_diagram_states[key]

    def __delitem__(self, key: int):
        del self._element_diagram_states[key]

    def __contains__(self, key: int):
        return key in self._element_diagram_states

    def generate_unnamed(self) -> int:
        \"\"\"
        Generate a number used in the name of an otherwise unnamed diagram
        \"\"\"
        self.unnamed_index += 1
        return self.unnamed_index

    def generate_index(self) -> int:
        \"\"\"
        Generate a number used to index a diagram
        \"\"\"
        self.index += 1
        return self.index

    def extract_into_diagram(self, el_id: int):
        \"\"\"
        Used when we encounter the same token twice in the same tree. When this
        happens, we replace all instances of that token with a terminal, and
        create a new subdiagram for the token
        \"\"\"
        position = self[el_id]

        # Replace the original definition of this element with a regular block
        if position.parent:
            ret = EditablePartial.from_call(railroad.NonTerminal, text=position.name)
            if \"item\" in position.parent.kwargs:
                position.parent.kwargs[\"item\"] = ret
            elif \"items\" in position.parent.kwargs:
                position.parent.kwargs[\"items\"][position.parent_index] = ret

        # If the element we're extracting is a group, skip to its content but keep the title
        if position.converted.func == railroad.Group:
            content = position.converted.kwargs[\"item\"]
        else:
            content = position.converted

        self.diagrams[el_id] = EditablePartial.from_call(
            NamedDiagram,
            name=position.name,
            diagram=EditablePartial.from_call(
                railroad.Diagram, content, **self.diagram_kwargs
            ),
            index=position.number,
        )

        del self[el_id]


def _worth_extracting(element: pyparsing.ParserElement) -> bool:
    \"\"\"
    Returns true if this element is worth having its own sub-diagram. Simply, if any of its children
    themselves have children, then its complex enough to extract
    \"\"\"
    children = element.recurse()
    return any(child.recurse() for child in children)


def _apply_diagram_item_enhancements(fn):
    \"\"\"
    decorator to ensure enhancements to a diagram item (such as results name annotations)
    get applied on return from _to_diagram_element (we do this since there are several
    returns in _to_diagram_element)
    \"\"\"

    def _inner(
        element: pyparsing.ParserElement,
        parent: typing.Optional[EditablePartial],
        lookup: ConverterState = None,
        vertical: int = None,
        index: int = 0,
        name_hint: str = None,
        show_results_names: bool = False,
        show_groups: bool = False,
    ) -> typing.Optional[EditablePartial]:

        ret = fn(
            element,
            parent,
            lookup,
            vertical,
            index,
            name_hint,
            show_results_names,
            show_groups,
        )

        # apply annotation for results name, if present
        if show_results_names and ret is not None:
            element_results_name = element.resultsName
            if element_results_name:
                # add \"*\" to indicate if this is a \"list all results\" name
                element_results_name += \"\" if element.modalResults else \"*\"
                ret = EditablePartial.from_call(
                    railroad.Group, item=ret, label=element_results_name
                )

        return ret

    return _inner


def _visible_exprs(exprs: Iterable[pyparsing.ParserElement]):
    non_diagramming_exprs = (
        pyparsing.ParseElementEnhance,
        pyparsing.PositionToken,
        pyparsing.And._ErrorStop,
    )
    return [
        e
        for e in exprs
        if not (e.customName or e.resultsName or isinstance(e, non_diagramming_exprs))
    ]


@_apply_diagram_item_enhancements
def _to_diagram_element(
    element: pyparsing.ParserElement,
    parent: typing.Optional[EditablePartial],
    lookup: ConverterState = None,
    vertical: int = None,
    index: int = 0,
    name_hint: str = None,
    show_results_names: bool = False,
    show_groups: bool = False,
) -> typing.Optional[EditablePartial]:
    \"\"\"
    Recursively converts a PyParsing Element to a railroad Element
    :param lookup: The shared converter state that keeps track of useful things
    :param index: The index of this element within the parent
    :param parent: The parent of this element in the output tree
    :param vertical: Controls at what point we make a list of elements vertical. If this is an integer (the default),
    it sets the threshold of the number of items before we go vertical. If True, always go vertical, if False, never
    do so
    :param name_hint: If provided, this will override the generated name
    :param show_results_names: bool flag indicating whether to add annotations for results names
    :returns: The converted version of the input element, but as a Partial that hasn't yet been constructed
    :param show_groups: bool flag indicating whether to show groups using bounding box
    \"\"\"
    exprs = element.recurse()
    name = name_hint or element.customName or element.__class__.__name__

    # Python's id() is used to provide a unique identifier for elements
    el_id = id(element)

    element_results_name = element.resultsName

    # Here we basically bypass processing certain wrapper elements if they contribute nothing to the diagram
    if not element.customName:
        if isinstance(
            element,
            (
                # pyparsing.TokenConverter,
                # pyparsing.Forward,
                pyparsing.Located,
            ),
        ):
            # However, if this element has a useful custom name, and its child does not, we can pass it on to the child
            if exprs:
                if not exprs[0].customName:
                    propagated_name = name
                else:
                    propagated_name = None

                return _to_diagram_element(
                    element.expr,
                    parent=parent,
                    lookup=lookup,
                    vertical=vertical,
                    index=index,
                    name_hint=propagated_name,
                    show_results_names=show_results_names,
                    show_groups=show_groups,
                )

    # If the element isn't worth extracting, we always treat it as the first time we say it
    if _worth_extracting(element):
        if el_id in lookup:
            # If we've seen this element exactly once before, we are only just now finding out that it's a duplicate,
            # so we have to extract it into a new diagram.
            looked_up = lookup[el_id]
            looked_up.mark_for_extraction(el_id, lookup, name=name_hint)
            ret = EditablePartial.from_call(railroad.NonTerminal, text=looked_up.name)
            return ret

        elif el_id in lookup.diagrams:
            # If we have seen the element at least twice before, and have already extracted it into a subdiagram, we
            # just put in a marker element that refers to the sub-diagram
            ret = EditablePartial.from_call(
                railroad.NonTerminal, text=lookup.diagrams[el_id].kwargs[\"name\"]
            )
            return ret

    # Recursively convert child elements
    # Here we find the most relevant Railroad element for matching pyparsing Element
    # We use ``items=[]`` here to hold the place for where the child elements will go once created
    if isinstance(element, pyparsing.And):
        # detect And's created with ``expr*N`` notation - for these use a OneOrMore with a repeat
        # (all will have the same name, and resultsName)
        if not exprs:
            return None
        if len(set((e.name, e.resultsName) for e in exprs)) == 1:
            ret = EditablePartial.from_call(
                railroad.OneOrMore, item=\"\", repeat=str(len(exprs))
            )
        elif _should_vertical(vertical, exprs):
            ret = EditablePartial.from_call(railroad.Stack, items=[])
        else:
            ret = EditablePartial.from_call(railroad.Sequence, items=[])
    elif isinstance(element, (pyparsing.Or, pyparsing.MatchFirst)):
        if not exprs:
            return None
        if _should_vertical(vertical, exprs):
            ret = EditablePartial.from_call(railroad.Choice, 0, items=[])
        else:
            ret = EditablePartial.from_call(railroad.HorizontalChoice, items=[])
    elif isinstance(element, pyparsing.Each):
        if not exprs:
            return None
        ret = EditablePartial.from_call(EachItem, items=[])
    elif isinstance(element, pyparsing.NotAny):
        ret = EditablePartial.from_call(AnnotatedItem, label=\"NOT\", item=\"\")
    elif isinstance(element, pyparsing.FollowedBy):
        ret = EditablePartial.from_call(AnnotatedItem, label=\"LOOKAHEAD\", item=\"\")
    elif isinstance(element, pyparsing.PrecededBy):
        ret = EditablePartial.from_call(AnnotatedItem, label=\"LOOKBEHIND\", item=\"\")
    elif isinstance(element, pyparsing.Group):
        if show_groups:
            ret = EditablePartial.from_call(AnnotatedItem, label=\"\", item=\"\")
        else:
            ret = EditablePartial.from_call(railroad.Group, label=\"\", item=\"\")
    elif isinstance(element, pyparsing.TokenConverter):
        ret = EditablePartial.from_call(
            AnnotatedItem, label=type(element).__name__.lower(), item=\"\"
        )
    elif isinstance(element, pyparsing.Opt):
        ret = EditablePartial.from_call(railroad.Optional, item=\"\")
    elif isinstance(element, pyparsing.OneOrMore):
        ret = EditablePartial.from_call(railroad.OneOrMore, item=\"\")
    elif isinstance(element, pyparsing.ZeroOrMore):
        ret = EditablePartial.from_call(railroad.ZeroOrMore, item=\"\")
    elif isinstance(element, pyparsing.Group):
        ret = EditablePartial.from_call(
            railroad.Group, item=None, label=element_results_name
        )
    elif isinstance(element, pyparsing.Empty) and not element.customName:
        # Skip unnamed \"Empty\" elements
        ret = None
    elif len(exprs) > 1:
        ret = EditablePartial.from_call(railroad.Sequence, items=[])
    elif len(exprs) > 0 and not element_results_name:
        ret = EditablePartial.from_call(railroad.Group, item=\"\", label=name)
    else:
        terminal = EditablePartial.from_call(railroad.Terminal, element.defaultName)
        ret = terminal

    if ret is None:
        return

    # Indicate this element's position in the tree so we can extract it if necessary
    lookup[el_id] = ElementState(
        element=element,
        converted=ret,
        parent=parent,
        parent_index=index,
        number=lookup.generate_index(),
    )
    if element.customName:
        lookup[el_id].mark_for_extraction(el_id, lookup, element.customName)

    i = 0
    for expr in exprs:
        # Add a placeholder index in case we have to extract the child before we even add it to the parent
        if \"items\" in ret.kwargs:
            ret.kwargs[\"items\"].insert(i, None)

        item = _to_diagram_element(
            expr,
            parent=ret,
            lookup=lookup,
            vertical=vertical,
            index=i,
            show_results_names=show_results_names,
            show_groups=show_groups,
        )

        # Some elements don't need to be shown in the diagram
        if item is not None:
            if \"item\" in ret.kwargs:
                ret.kwargs[\"item\"] = item
            elif \"items\" in ret.kwargs:
                # If we've already extracted the child, don't touch this index, since it's occupied by a nonterminal
                ret.kwargs[\"items\"][i] = item
                i += 1
        elif \"items\" in ret.kwargs:
            # If we're supposed to skip this element, remove it from the parent
            del ret.kwargs[\"items\"][i]

    # If all this items children are none, skip this item
    if ret and (
        (\"items\" in ret.kwargs and len(ret.kwargs[\"items\"]) == 0)
        or (\"item\" in ret.kwargs and ret.kwargs[\"item\"] is None)
    ):
        ret = EditablePartial.from_call(railroad.Terminal, name)

    # Mark this element as \"complete\", ie it has all of its children
    if el_id in lookup:
        lookup[el_id].complete = True

    if el_id in lookup and lookup[el_id].extract and lookup[el_id].complete:
        lookup.extract_into_diagram(el_id)
        if ret is not None:
            ret = EditablePartial.from_call(
                railroad.NonTerminal, text=lookup.diagrams[el_id].kwargs[\"name\"]
            )

    return ret

"""

import os
import types
import zipfile
import sys
import io

class ZipImporter(object):
    def __init__(self, zip_file):
        self.zfile = zip_file
        self._paths = [x.filename for x in self.zfile.filelist]
        
    def _mod_to_paths(self, fullname):
        # get the python module name
        py_filename = fullname.replace(".", os.sep) + ".py"
        # get the filename if it is a package/subpackage
        py_package = fullname.replace(".", os.sep) + "/__init__.py"
        if py_filename in self._paths:
            return py_filename
        elif py_package in self._paths:
            return py_package
        else:
            return None

    def find_module(self, fullname, path):
        if self._mod_to_paths(fullname) is not None:
            return self
        return None

    def load_module(self, fullname):
        filename = self._mod_to_paths(fullname)
        if not filename in self._paths:
            raise ImportError(fullname)
        new_module = types.ModuleType(fullname)
        sys.modules[fullname]=new_module
        if filename.endswith("__init__.py"):
            new_module.__path__ = [] 
            new_module.__package__ = fullname
        else:
            new_module.__package__ = fullname.rpartition('.')[0]
        exec(self.zfile.open(filename, 'r').read(),new_module.__dict__)
        new_module.__file__ = filename
        new_module.__loader__ = self
        return new_module

module_zip=zipfile.ZipFile(io.BytesIO(),"w")
for key in module_dict:
    module_zip.writestr(key,module_dict[key])

module_importer=ZipImporter(module_zip)
sys.meta_path.insert(0,module_importer)

from pyparsing import *

sys.meta_path.remove(module_importer)

for key in sys.modules.copy():
    if key=="pyparsing" or key.startswith("pyparsing."):
        del sys.modules[key]
