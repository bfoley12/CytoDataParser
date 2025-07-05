import re
import operator
from typing import Any, Callable


def greater_than(threshold: float) -> Callable[[Any], bool]:
    return lambda x: x is not None and x > threshold


def less_than(threshold: float) -> Callable[[Any], bool]:
    return lambda x: x is not None and x < threshold


def in_list(values: list) -> Callable[[Any], bool]:
    return lambda x: x in values


def not_null() -> Callable[[Any], bool]:
    return lambda x: x is not None


def matches_regex(pattern: str) -> Callable[[Any], bool]:
    compiled = re.compile(pattern)
    return lambda x: isinstance(x, str) and compiled.search(x) is not None

# TODO: allow for str == int, ie MouseID == 1
# TODO: fix return types: str is returned in parse_single, which is not callable


def parse_string_condition(expr: str) -> Callable[[Any], bool]:
    """
    Convert string expressions like '> 10', '<= 50', or compound forms like '> 10 and < 20'
    into a callable predicate.
    """
    ops = {
        '>=': operator.ge,
        '<=': operator.le,
        '!=': operator.ne,
        '==': operator.eq,
        '>': operator.gt,
        '<': operator.lt
    }

    def parse_single(expression: str) -> Callable[[Any], bool]:
        for op_str in sorted(ops.keys(), key=len, reverse=True):
            if expression.strip().startswith(op_str):
                value_str = expression.strip()[len(op_str):].strip()
                try:
                    value = eval(value_str)
                except Exception:
                    value = value_str
                return lambda x: ops[op_str](x, value)
        # If no known operator matched, return a predicate that always fails
        return lambda x: False
        
    if ' and ' in expr:
        parts = expr.split(' and ')
        funcs = [parse_single(part.strip()) for part in parts]
        return lambda x: all(f(x) for f in funcs)

    return parse_single(expr.strip())



def from_range(r: range) -> Callable[[Any], bool]:
    """
    Convert a Python range object to a predicate.
    """
    return lambda x: x in r
