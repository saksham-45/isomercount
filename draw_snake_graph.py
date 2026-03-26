#!/usr/bin/env python3
"""Compatibility wrapper.

Implementation: src/snakes/draw_snake_graph.py

Note: tests import underscore-prefixed helpers from this module, so we re-export
the implementation module namespace instead of using ``from ... import *`` (which
would hide underscore names).
"""

from src.snakes import draw_snake_graph as _impl  # type: ignore

# Re-export everything except dunder names.
for _k, _v in _impl.__dict__.items():
    if not _k.startswith("__"):
        globals()[_k] = _v


if __name__ == "__main__":
    _impl.main()
