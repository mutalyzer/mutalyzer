from importlib.metadata import metadata


def _get_metadata(name):
    meta = metadata(__package__)
    return meta.get(name, "")


_copyright_notice = "Copyright (c) {} <{}>".format(
    _get_metadata("Author"), _get_metadata("Author-email")
)

usage = [_get_metadata("Summary"), _copyright_notice]


def doc_split(func):
    return func.__doc__.split("\n\n")[0]


def version(name):
    return "{} version {}\n\n{}\nHomepage: {}".format(
        _get_metadata("Name"),
        _get_metadata("Version"),
        _copyright_notice,
        _get_metadata("Home-page"),
    )
