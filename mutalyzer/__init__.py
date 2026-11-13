from email.utils import parseaddr
from importlib.metadata import PackageNotFoundError, metadata
from re import split


def _extract(key: str, delim: str = r"[^\s\S]", index: int = 0) -> str:
    try:
        value = metadata(__package__).get(key, "")
    except PackageNotFoundError:
        return "<NO DATA>"
    return split(delim, value)[index]


def doc_split(func):
    return func.__doc__.split("\n\n")[0]


_url = _extract("Project-URL", r", ", 1)
_year = "2019-2026"

project = _extract("Name")
release = _extract("Version")
author, _email = parseaddr(_extract("Author-email"))
# Not `description`, which is the mutalyzer.description module once imported.
summary = _extract("Summary")
copyright = f"{_year} by {author} <{_email}>"
info = f"""{project} version {release}

Copyright (c) {copyright}
Homepage: {_url}"""
