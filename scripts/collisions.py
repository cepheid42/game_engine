from dataclasses import dataclass

@dataclass
class Collision:
    groups: tuple = ()
    types: tuple = ()
    coulomb_log: float = -1.0
    rate_mult: float = 1.0
    self_scatter: bool = False
    step_interval: int = 1

    def __repr__(self, indent=0):
        return (f'{indent * ' '}std::tuple('
                f'"{self.groups[0]}", "{self.groups[1]}", '
                f'{self.coulomb_log}, '
                f'{self.rate_mult}, '
                f'{self.step_interval}, '
                f'{str(self.self_scatter).lower()}),')

