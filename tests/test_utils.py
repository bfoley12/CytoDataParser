from cytodataparser.utils.predicates import parse_string_condition, from_range, matches_regex

def test_string_condition_parsing():
    cond = parse_string_condition("> 10 and < 20")
    assert cond(15)
    assert not cond(25)

    cond = parse_string_condition("== 10")
    assert cond(10)
    assert not cond(15)

    cond = parse_string_condition("<= 10")
    assert cond(9)
    assert not cond(15)

    cond = parse_string_condition(">= 10")
    assert not cond(9)
    assert cond(15)

    cond = parse_string_condition("!= 10")
    assert cond(9)
    assert not cond(10)

    cond = parse_string_condition("!= Hello")
    assert not cond("Hello")
    assert cond("Goodbye")

def test_range_condition():
    r = from_range(range(5, 10))
    assert r(6)
    assert not r(4)

def test_regex_condition():
    cond = matches_regex(r"^AB\d{3}Z$")
    assert cond("AB123Z")
    assert not cond("AB12Z")