#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Time-domain induced polarization (TDIP) Data Manager."""

from .tdip import TDIP
from .hirip import HIRIP
from .decay import Decay
from .modelling import (DCIPMModelling, DCIPSeigelModelling,
                        ColeColeTD, DCIPMSmoothModelling,
                        CCTDModelling, MultiDebyeTDModelling)

TDIPdata = TDIP  # backward compatibility

__all__ = ['TDIP', 'TDIPdata', 'HIRIP', 'Decay', 'DCIPMModelling',
           'DCIPSeigelModelling', 'ColeColeTD', 'MultiDebyeTDModelling']
