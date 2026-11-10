#!/usr/bin/env python3
"""Compatibility entry point; new development uses sbml2map.py."""
from sbml2map import BacterialStrategy, discovery_context, main, parser, run

__all__ = ['BacterialStrategy', 'discovery_context', 'main', 'parser', 'run']

if __name__ == '__main__':
    main()
