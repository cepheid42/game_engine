from dataclasses import dataclass

import numpy as np

from .header_utils import *

@dataclass
class TimeParams:
    dt : float
    nt : int
    t_end : float

    def __init__(self, dt : float, nt : None | int = None, t_end : None | float = None):
        self.dt = dt
        if nt is None and t_end is not None:
            self.nt = int(np.ceil(t_end / dt))
            self.t_end = t_end
        elif nt is not None and t_end is None:
            self.nt = nt
            self.t_end = dt * nt
        else:
            raise(RuntimeError)


    def __repr__(self) -> str:
        return (
            section_label('Time') +
            parameter_string('auto', 'dt', self.dt) +
            parameter_string('auto', 'Nt', str(self.nt) + 'lu') +
            parameter_string('auto', 't_end', self.t_end)
        )
