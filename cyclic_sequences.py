#!/usr/bin/env python3
"""Compatibility wrapper.

The implementation lives in src/cyclic/cyclic_sequences.py.
This wrapper keeps existing imports/CLI invocations working:
  python3 cyclic_sequences.py --max-n 20
  from cyclic_sequences import count_adjacency_burnside
"""

from src.cyclic.cyclic_sequences import *  # noqa: F401,F403
from src.cyclic.cyclic_sequences import main


if __name__ == "__main__":
    main()
